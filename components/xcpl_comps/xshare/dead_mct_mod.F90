module dead_mct_mod

  use esmf            , only : esmf_clock
  use shr_kind_mod    , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : shr_const_pi
  use mct_mod         , only : mct_gsmap, mct_ggrid, mct_avect, mct_ggrid_init, mct_gsmap_lsize, mct_ggrid_lsize
  use mct_mod         , only : mct_avect_lsize, MCT_AVECT_NRATTR, mct_avect_indexra,mct_avect_zero, mct_aVect_exportRList2c
  use mct_mod         , only : mct_ggrid_importiattr, mct_ggrid_importrattr, mct_gsmap_init, mct_aVect_init
  use mct_mod         , only : mct_gsmap_orderedpoints
  use dead_data_mod   , only : dead_grid_lat, dead_grid_lon, dead_grid_area, dead_grid_mask, dead_grid_frac, dead_grid_index
  use dead_mod        , only : dead_setnewgrid, dead_read_inparms
  use seq_flds_mod    , only : seq_flds_dom_coord, seq_flds_dom_other
  use seq_timemgr_mod , only : seq_timemgr_EClockGetData
  use iso_c_binding

  implicit none
  private
  save

  public  :: dead_init_mct, dead_run_mct, dead_final_mct
  private :: dead_domain_mct
#ifdef HAVE_MOAB
  public :: dead_init_moab
  private :: define_reset_fields_moab
  ! Saved grid dimensions per model index (1=atm..7=wav)
  ! Used by dead_run_mct to reorder data for MOAB entity ordering
  integer(IN) :: moab_nxg(7) = 0, moab_nyg(7) = 0
#endif

!===============================================================================
contains
!===============================================================================

  subroutine dead_init_mct(model, Eclock, x2d, d2x, &
         flds_x2d, flds_d2x, &
         gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, nxg, nyg, flood)

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in)    :: model
    type(ESMF_Clock) , intent(inout) :: EClock
    type(mct_aVect)  , intent(inout) :: x2d         ! driver -> dead
    type(mct_aVect)  , intent(inout) :: d2x         ! dead   -> driver
    character(len=*) , intent(in)    :: flds_x2d
    character(len=*) , intent(in)    :: flds_d2x
    type(mct_gsMap)  , pointer       :: gsMap       ! model global sep map (output)
    type(mct_gGrid)  , pointer       :: ggrid       ! model ggrid (output)
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid (output)
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)      , intent(in)    :: compid      ! mct comp id
    integer(IN)      , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in)    :: master_task ! task number of master task
    integer(IN)      , intent(in)    :: inst_index  ! instance number
    character(len=*) , intent(in)    :: inst_suffix ! char string associated with instance
    character(len=*) , intent(in)    :: inst_name   ! fullname of current instance (ie. "lnd_0001")
    integer(IN)      , intent(in)    :: logunit     ! logging unit number
    integer(IN)      , intent(out)   :: nxg         ! global dim i-direction
    integer(IN)      , intent(out)   :: nyg         ! global dim j-direction
    logical, intent(out),optional   :: flood

    !--- local variables ---
    integer(IN)              :: ierr          ! error code
    integer(IN)              :: local_comm    ! local communicator
    integer(IN)              :: mype          ! pe info
    integer(IN)              :: totpe         ! total number of pes
    integer(IN), allocatable :: gindex(:)     ! global index
    integer(IN)              :: lsize
    integer(IN)              :: nproc_x
    integer(IN)              :: seg_len
    integer(IN)              :: decomp_type
    logical                  :: lflood=.false.

    !--- formats ---
    character(*), parameter :: F00   = "('(',a,'_init_mct) ',8a)"
    character(*), parameter :: F01   = "('(',a,'_init_mct) ',a,4i8)"
    character(*), parameter :: F02   = "('(',a,'_init_mct) ',a,4es13.6)"
    character(*), parameter :: F90   = "('(',a,'_init_mct) ',73('='))"
    character(*), parameter :: F91   = "('(',a,'_init_mct) ',73('-'))"
    character(*), parameter :: subName = "(dead_init_mct) "
    !-------------------------------------------------------------------------------

    ! Determine communicator groups and sizes

    local_comm = mpicom
    call MPI_COMM_RANK(local_comm,mype ,ierr)
    call MPI_COMM_SIZE(local_comm,totpe,ierr)

    ! Read input parms

    call dead_read_inparms(model, mpicom, my_task, master_task, &
       inst_index, inst_suffix, inst_name, logunit, &
       nxg, nyg, decomp_type, nproc_x, seg_len, lflood)

    if (present(flood)) then
      flood=lflood
    endif

    ! Initialize grid

    call dead_setNewGrid(decomp_type, nxg, nyg, totpe, mype, logunit, &
                         lsize, gbuf, seg_len, nproc_x)

    ! Initialize MCT global seg map

    allocate(gindex(lsize))
    gindex(:) = gbuf(:,dead_grid_index)
    call mct_gsMap_init( gsMap, gindex, mpicom, compid, lsize, nxg*nyg )
    deallocate(gindex)

    ! Initialize MCT domain

    call dead_domain_mct(mpicom, gbuf, gsMap, logunit, ggrid)

    ! Initialize MCT attribute vectors

    call mct_aVect_init(d2x, rList=flds_d2x, lsize=lsize)
    call mct_aVect_zero(d2x)

    call mct_aVect_init(x2d, rList=flds_x2d, lsize=lsize)
    call mct_aVect_zero(x2d)

  end subroutine dead_init_mct

  !===============================================================================
  subroutine dead_run_mct(model, EClock, x2d, d2x, &
       gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, logunit, mbdomain)

#ifdef HAVE_MOAB
      use iMOAB      , only : iMOAB_SetDoubleTagStorage
#endif
    implicit none

    ! !DESCRIPTION: run method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in)    :: model
    type(ESMF_Clock) , intent(inout) :: EClock
    type(mct_aVect)  , intent(inout) :: x2d         ! driver -> dead
    type(mct_aVect)  , intent(inout) :: d2x         ! dead   -> driver
    type(mct_gsMap)  , pointer       :: gsMap       ! model global sep map (output)
    type(mct_gGrid)  , pointer       :: ggrid       ! model ggrid (output)
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)      , intent(in)    :: compid      ! mct comp id
    integer(IN)      , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in)    :: master_task ! task number of master task
    integer(IN)      , intent(in)    :: logunit     ! logging unit number
    integer(IN)      , intent(in), optional    :: mbdomain    ! MOAB domain handle (only useful for MOAB coupler)
#ifdef HAVE_MOAB
    integer(IN)      :: ierr        ! error code
    integer(IN)      :: ent_type    ! entity type: 0=vertex, 1=element
    integer(IN)      :: ntri_run, tri_cnt, quad_cnt, j_run, gid_run
    real(R8), allocatable :: data_ordered(:,:)  ! reordered data for MOAB
#endif

    !--- local ---
    integer(IN)             :: CurrentYMD        ! model date
    integer(IN)             :: CurrentTOD        ! model sec into model date
    integer(IN)             :: n                 ! index
    integer(IN)             :: nf                ! fields loop index
    integer(IN)             :: ki                ! index of ifrac
    integer(IN)             :: lsize             ! size of AttrVect
    real(R8)                :: lat               ! latitude
    real(R8)                :: lon               ! longitude

    real(R8), pointer    :: data_d2x(:,:) ! data_x2d(:,:)
    integer                 :: nflds_d2x, nflds_x2d
    character(len=2048)    :: flds_x2d
    character(len=2048)    :: flds_d2x
    integer                 :: ncomp
    character(*), parameter :: F04   = "('(',a,'_run_mct) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dead_run_mct) "
    !-------------------------------------------------------------------------------

    ! PACK (currently no unpacking)

    selectcase(model)
    case('atm')
       ncomp = 1
    case('lnd')
       ncomp = 2
    case('ice')
       ncomp = 3
    case('ocn')
       ncomp = 4
    case('glc')
       ncomp = 5
    case('rof')
       ncomp = 6
    case('wav')
       ncomp = 7
    end select

#ifdef HAVE_MOAB
    ! Set entity type for MOAB tag operations: elements for all dead components
    ent_type = 1  ! elements (all dead components have full RLL mesh)
#endif

    lsize = mct_avect_lsize(x2d)
    nflds_d2x = mct_avect_nRattr(d2x)
    nflds_x2d = mct_avect_nRattr(x2d)

#ifdef HAVE_MOAB
    flds_d2x = trim(mct_aVect_exportRList2c(d2x)) ! exportRListToChar
    !flds_x2d = trim(mct_aVect_exportRList2c(x2d)) ! exportRListToChar

    allocate(data_d2x(lsize,nflds_d2x))
   !  allocate(data_x2d(lsize,nflds_x2d))
#endif

    if (model.eq.'rof') then

       do nf=1,nflds_d2x
          do n=1,lsize
             d2x%rAttr(nf,n) = (nf+1) * 1.0_r8
#ifdef HAVE_MOAB
             data_d2x(n,nf) = d2x%rAttr(nf,n)
#endif
          enddo
       enddo

    else if (model.eq.'glc') then

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x%rAttr(nf,n) = (nf*100)                   &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin (SHR_CONST_PI*lon/180.0_R8)       &
                  *  sin (SHR_CONST_PI*lon/180.0_R8)       &
                  + (ncomp*10.0_R8)
#ifdef HAVE_MOAB
               data_d2x(n,nf) = d2x%rAttr(nf,n)
#endif
          enddo
       enddo

    else

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x%rAttr(nf,n) = (nf*100)                   &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin((SHR_CONST_PI*lon/180.0_R8)       &
                  -      (ncomp-1)*(SHR_CONST_PI/3.0_R8) ) &
                  + (ncomp*10.0_R8)
#ifdef HAVE_MOAB
               data_d2x(n,nf) = d2x%rAttr(nf,n)
#endif
          enddo
       enddo

    endif

    selectcase(model)
    case('ice')

       ki = mct_aVect_indexRA(d2x,"Si_ifrac",perrWith=subname)
       d2x%rAttr(ki,:) = min(1.0_R8,max(0.0_R8,d2x%rAttr(ki,:)))
#ifdef HAVE_MOAB
       data_d2x(:,ki) = d2x%rAttr(ki,:)
#endif

    case('glc')

       ki = mct_aVect_indexRA(d2x,"Sg_icemask",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8
#ifdef HAVE_MOAB
       data_d2x(:,ki) = d2x%rAttr(ki,:)
#endif

       ki = mct_aVect_indexRA(d2x,"Sg_icemask_coupled_fluxes",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8
#ifdef HAVE_MOAB
       data_d2x(:,ki) = d2x%rAttr(ki,:)
#endif

       ki = mct_aVect_indexRA(d2x,"Sg_ice_covered",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8
#ifdef HAVE_MOAB
       data_d2x(:,ki) = d2x%rAttr(ki,:)
#endif
    end select

#ifdef HAVE_MOAB
    ! Import data into MOAB structures
    ! For fullmesh (non-ATM), reorder data to match MOAB entity iteration order
    ! (all MBTRI first, then all MBQUAD)
    if (ent_type == 1 .and. moab_nxg(ncomp) > 0) then
       ! Count triangle elements on this partition
       ntri_run = 0
       do n = 1, lsize
          gid_run = int(gbuf(n, dead_grid_index))
          j_run = (gid_run - 1) / moab_nxg(ncomp)
          if (j_run == 0 .or. j_run == moab_nyg(ncomp) - 1) &
             ntri_run = ntri_run + 1
       enddo
       ! Reorder: tris first, then quads
       allocate(data_ordered(lsize, nflds_d2x))
       tri_cnt = 0
       quad_cnt = 0
       do n = 1, lsize
          gid_run = int(gbuf(n, dead_grid_index))
          j_run = (gid_run - 1) / moab_nxg(ncomp)
          if (j_run == 0 .or. j_run == moab_nyg(ncomp) - 1) then
             tri_cnt = tri_cnt + 1
             data_ordered(tri_cnt, :) = data_d2x(n, :)
          else
             quad_cnt = quad_cnt + 1
             data_ordered(ntri_run + quad_cnt, :) = data_d2x(n, :)
          endif
       enddo
       ierr = iMOAB_SetDoubleTagStorage ( mbdomain, trim(flds_d2x) // C_NULL_CHAR, lsize*nflds_d2x, &
                                        ent_type, &
                                        data_ordered)
       deallocate(data_ordered)
    else
       ierr = iMOAB_SetDoubleTagStorage ( mbdomain, trim(flds_d2x) // C_NULL_CHAR, lsize*nflds_d2x, &
                                        ent_type, &
                                        data_d2x)
    endif
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set flds_d2x fields tag ')

    deallocate(data_d2x)
    !deallocate(data_x2d)
#endif

    ! log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
       write(logunit,F04) model,trim(model),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if

  end subroutine dead_run_mct

  !===============================================================================
  subroutine dead_final_mct(model, my_task, master_task, logunit)

    ! !DESCRIPTION: finalize method for datm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in) :: model
    integer(IN)      , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in) :: master_task ! task number of master task
    integer(IN)      , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dead_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dead_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dead_comp_final) "
    !-------------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(model),': end of main integration loop'
       write(logunit,F91)
    end if

  end subroutine dead_final_mct

  !===============================================================================
  subroutine dead_domain_mct( mpicom, gbuf, gsMap, logunit, domain )

    !-------------------------------------------------------------------
    !---arguments---
    integer(IN)    , intent(in)  :: mpicom
    real(R8)       , intent(in)  :: gbuf(:,:)
    type(mct_gsMap), intent(in)  :: gsMap
    integer(IN)    , intent(in)  :: logunit
    type(mct_ggrid), intent(out) :: domain

    !---local variables---
    integer(IN)          :: my_task     ! mpi task within communicator
    integer(IN)          :: lsize       ! temporary
    integer(IN)          :: n    ,j,i   ! indices
    integer(IN)          :: ier         ! error status
    real(R8), pointer    :: data(:)     ! temporary
    integer(IN), pointer :: idata(:)    ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct dead domain
    !
    call mct_gGrid_init( GGrid=domain, CoordChars=trim(seq_flds_dom_coord), &
         OtherChars=trim(seq_flds_dom_other), &
         lsize=mct_gsMap_lsize(gsMap, mpicom) )
    call mct_aVect_zero(domain%data)
    !
    ! Allocate memory
    !
    lsize = mct_gGrid_lsize(domain)
    if (size(gbuf,dim=1) /= lsize) then
       call shr_sys_abort('mct_dead_domain size error')
    endif
    allocate(data(lsize))
    allocate(idata(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mpi_comm_rank(mpicom, my_task, ier)
    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
    call mct_gGrid_importIAttr(domain,'GlobGridNum',idata,lsize)
    !
    call mct_aVect_zero(domain%data)
    data(:) = -9999.0_R8  ! generic special value
    call mct_gGrid_importRAttr(domain,"lat" ,data,lsize)
    call mct_gGrid_importRAttr(domain,"lon" ,data,lsize)
    call mct_gGrid_importRAttr(domain,"area",data,lsize)
    call mct_gGrid_importRAttr(domain,"frac",data,lsize)

    data(:) = 0.0_R8  ! generic special value
    call mct_gGrid_importRAttr(domain,"mask" ,data,lsize)
    call mct_gGrid_importRAttr(domain,"aream",data,lsize)
    !
    ! Fill in correct values for domain components
    !
    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_lat)
    enddo
    call mct_gGrid_importRAttr(domain,"lat",data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_lon)
    enddo
    call mct_gGrid_importRAttr(domain,"lon",data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_area)
    enddo
    call mct_gGrid_importRAttr(domain,"area",data,lsize)
    call mct_gGrid_importRAttr(domain,"aream",data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_mask)
    enddo
    call mct_gGrid_importRAttr(domain,"mask"   ,data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_frac)
    enddo
    call mct_gGrid_importRAttr(domain,"frac"   ,data,lsize)

    deallocate(data)
    deallocate(idata)

  end subroutine dead_domain_mct
  !===============================================================================

#ifdef HAVE_MOAB
  subroutine dead_init_moab( mbdomain, model, gsMap, gbuf, x2d, d2x, mpicom, compid, logunit, nxg, nyg )

    use iMOAB, only: iMOAB_CreateVertices, &
                     iMOAB_CreateElements, &
                     iMOAB_DefineTagStorage, &
                     iMOAB_SetIntTagStorage, &
                     iMOAB_SetDoubleTagStorage, &
                     iMOAB_ResolveSharedEntities, &
                     iMOAB_UpdateMeshInfo, &
                     iMOAB_WriteMesh
      implicit none
    !-------------------------------------------------------------------
    !---arguments---
    integer(IN)    , intent(inout)   :: mbdomain
    character(len=*) , intent(in) :: model
    type(mct_gsMap), intent(in)   :: gsMap
    real(R8)       , intent(in)   :: gbuf(:,:)
    type(mct_aVect)  , intent(in) :: x2d         ! driver -> dead
    type(mct_aVect)  , intent(in) :: d2x         ! dead   -> driver
    integer(IN)    , intent(in)   :: mpicom
    integer(IN)    , intent(in)   :: compid      ! mct comp id
    integer(IN)    , intent(in)   :: logunit
    integer(IN)    , intent(in)   :: nxg          ! global grid size in x-direction
    integer(IN)    , intent(in)   :: nyg          ! global grid size in y-direction

    !---local variables---
    integer(IN)          :: ierr        ! error code
    integer(IN)          :: my_task     ! mpi task within communicator
    integer(IN)          :: lsize       ! temporary
    integer(IN)          :: n    ,j,i   ! indices
    integer(IN)          :: ier         ! error status
    real(R8)             :: lonv, latv  ! lon/lat values
    character(CL)        :: tagname
    integer(IN)          :: tagindex, etype
    real(R8), pointer    :: data(:)     ! temporary
    integer(IN), pointer :: idata(:)    ! temporary
    real(R8), dimension(:), allocatable :: moab_vert_coords  ! temporary
    integer(IN), allocatable :: entity_order(:)  ! maps MOAB elem pos -> gbuf index
    integer(IN)          :: midx        ! model index for moab_nxg/nyg
    logical          :: fullmesh
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
    character(*), parameter :: subName = "(dead_init_moab) "
    !-------------------------------------------------------------------
    !
    ! Initialize MOAB dead domain
    !

    ! Allocate memory
    !
    lsize = size(gbuf,dim=1)
    !
    ! Initialize attribute vector with special value
    !
    call mpi_comm_rank(mpicom, my_task, ier)

    ! Remove the previously created vertex cloud; we will create a proper RLL mesh with vertices at corners
    ! This is done by deleting the mesh and starting fresh with corner-based vertices

    fullmesh = .true.
    etype = 1  ! elements (all dead components create full RLL mesh)

    ! Save grid dimensions per model for use in dead_run_mct
    selectcase(model)
    case('atm'); midx = 1
    case('lnd'); midx = 2
    case('ice'); midx = 3
    case('ocn'); midx = 4
    case('glc'); midx = 5
    case('rof'); midx = 6
    case('wav'); midx = 7
    case default; midx = 0
    end select
    if (midx > 0) then
       moab_nxg(midx) = nxg
       moab_nyg(midx) = nyg
    endif

    if (fullmesh) then
       ! Create RLL mesh with proper periodic longitude wrapping and polar vertices
       ! For nxg x nyg element grid on a sphere:
       !   - South pole: single vertex at lat=-90
       !   - North pole: single vertex at lat=+90
       !   - Interior vertices: nxg per latitude row (periodic in longitude)
       !   - Pole row elements (j=0, j=nyg-1): triangles
       !   - Interior elements: quads
       !
       ! This implementation scales with local element count (O(lsize * log(lsize)))
       ! and does not allocate or iterate over global-sized arrays.
       block
          integer(IN)          :: nverts_local   ! number of unique local vertices
          integer(IN)          :: nelems         ! number of elements
          integer(IN), allocatable :: vert_gid(:)      ! sorted unique vertex global IDs
          real(R8), allocatable    :: vert_coords(:)   ! vertex coords in 3D (x,y,z)
          integer(IN)          :: iv, ix, iy, vgid
          integer(IN)          :: elem_type
          real(R8)             :: dlat, dlon, lat_v, lon_v
          integer(IN)          :: i_elem, j_elem  ! global element indices (0-based)
          integer(IN)          :: elem_gid        ! global element ID (1-based)
          integer(IN)          :: ntri_elems, nquad_elems, tri_idx, quad_idx
          integer(IN), allocatable :: tri_conn(:), quad_conn(:)
          integer(IN)          :: ix_r            ! right-side ix with wrapping
          integer(IN)          :: vgid_sp, vgid_np ! south/north pole vertex global IDs
          integer(IN)          :: vgid_bl, vgid_br, vgid_tl, vgid_tr
          integer(IN), allocatable :: elem_gids(:)  ! element global IDs for MOAB tag
          integer(IN), allocatable :: all_vgids(:)  ! workspace for collecting vertex GIDs
          integer(IN)          :: nv_total, k      ! workspace counters

          ! Grid spacing for regular lat-lon mesh
          ! Note: gbuf longitudes start at 0 deg (from dead_setNewGrid)
          dlat = 180.0_R8 / real(nyg, R8)
          dlon = 360.0_R8 / real(nxg, R8)

          nelems = lsize

          ! Global vertex ID scheme for periodic RLL mesh:
          !   South pole vertex:  vgid_sp = 1
          !   Non-pole vertex at column ix (0..nxg-1), row iy (1..nyg-1):
          !     vgid = 2 + (iy - 1) * nxg + ix
          !   North pole vertex:  vgid_np = 2 + (nyg - 1) * nxg
          ! Total unique vertices = (nyg - 1) * nxg + 2
          vgid_sp = 1
          vgid_np = 2 + (nyg - 1) * nxg

          ! --- Step 1: Collect all vertex GIDs referenced by local elements ---
          ! Each element has at most 4 vertices, so 4*lsize is an upper bound.
          allocate(all_vgids(4 * lsize))
          nv_total = 0
          ntri_elems = 0
          nquad_elems = 0
          do n = 1, lsize
             elem_gid = int(gbuf(n, dead_grid_index))
             i_elem = mod(elem_gid - 1, nxg)      ! 0-based column
             j_elem = (elem_gid - 1) / nxg         ! 0-based row
             ix_r = mod(i_elem + 1, nxg)            ! right column with wrapping

             if (j_elem == 0) then
                ! South pole triangle: pole + 2 top-edge vertices at iy=1
                ntri_elems = ntri_elems + 1
                nv_total = nv_total + 1; all_vgids(nv_total) = vgid_sp
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + i_elem
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + ix_r
             else if (j_elem == nyg - 1) then
                ! North pole triangle: 2 bottom-edge vertices at iy=nyg-1 + pole
                ntri_elems = ntri_elems + 1
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + (nyg-2)*nxg + i_elem
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + (nyg-2)*nxg + ix_r
                nv_total = nv_total + 1; all_vgids(nv_total) = vgid_np
             else
                ! Interior quad: 4 corner vertices
                nquad_elems = nquad_elems + 1
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + (j_elem-1)*nxg + i_elem
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + (j_elem-1)*nxg + ix_r
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + j_elem*nxg + ix_r
                nv_total = nv_total + 1; all_vgids(nv_total) = 2 + j_elem*nxg + i_elem
             endif
          enddo

          ! --- Step 2: Sort and deduplicate to get unique local vertex GIDs ---
          ! Uses heapsort for guaranteed O(n log n) and no extra memory.
          call dead_heapsort_int(all_vgids, nv_total)

          ! Compact: remove duplicates (array is now sorted)
          nverts_local = 1
          do k = 2, nv_total
             if (all_vgids(k) /= all_vgids(nverts_local)) then
                nverts_local = nverts_local + 1
                all_vgids(nverts_local) = all_vgids(k)
             endif
          enddo

          ! Move unique GIDs into final array
          allocate(vert_gid(nverts_local))
          vert_gid(1:nverts_local) = all_vgids(1:nverts_local)
          deallocate(all_vgids)

          ! --- Step 3: Compute vertex coordinates only for local vertices ---
          ! Coordinates are derived analytically from the vertex global ID.
          allocate(vert_coords(nverts_local * 3))
          do iv = 1, nverts_local
             vgid = vert_gid(iv)
             if (vgid == vgid_sp) then
                ! South pole: lat = -90 deg
                lat_v = -0.5_R8 * shr_const_pi
                lon_v = 0.0_R8
             else if (vgid == vgid_np) then
                ! North pole: lat = +90 deg
                lat_v = 0.5_R8 * shr_const_pi
                lon_v = 0.0_R8
             else
                ! Non-pole vertex: decode ix, iy from vgid = 2 + (iy-1)*nxg + ix
                ix = mod(vgid - 2, nxg)
                iy = (vgid - 2) / nxg + 1
                ! Longitude starts at 0 deg (matching dead_setNewGrid convention)
                lat_v = (-90.0_R8 + real(iy, R8) * dlat) * shr_const_pi / 180.0_R8
                lon_v = real(ix, R8) * dlon * shr_const_pi / 180.0_R8
             endif
             vert_coords(3*iv-2) = cos(lat_v) * cos(lon_v)
             vert_coords(3*iv-1) = cos(lat_v) * sin(lon_v)
             vert_coords(3*iv  ) = sin(lat_v)
          enddo

          ! Create vertices in MOAB
          ierr = iMOAB_CreateVertices(mbdomain, nverts_local*3, 3, vert_coords)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to create vertices for RLL mesh')

          ! --- Step 4: Build element connectivity ---
          ! Uses binary search in vert_gid(:) instead of a global-sized lookup table.

          allocate(tri_conn(max(ntri_elems * 3, 1)))
          allocate(quad_conn(max(nquad_elems * 4, 1)))

          ! Build element connectivity using 1-based local vertex indices
          ! (iMOAB_CreateElements expects Fortran 1-based connectivity)
          ! Also build element GLOBAL_ID array ordered as: tris first, then quads
          ! (matching MOAB's entity iteration order: MBTRI before MBQUAD)
          ! entity_order maps MOAB element position -> local gbuf index
          allocate(elem_gids(nelems))
          allocate(entity_order(nelems))
          tri_idx = 0
          quad_idx = 0
          do n = 1, lsize
             elem_gid = int(gbuf(n, dead_grid_index))
             i_elem = mod(elem_gid - 1, nxg)
             j_elem = (elem_gid - 1) / nxg
             ix_r = mod(i_elem + 1, nxg)

             if (j_elem == 0) then
                ! South pole triangle: SP, (ix_r, iy=1), (i_elem, iy=1)
                ! Winding order gives outward normal toward -z (correct at south pole)
                tri_idx = tri_idx + 1
                tri_conn(3*(tri_idx-1) + 1) = dead_bsearch_int(vert_gid, nverts_local, vgid_sp)
                tri_conn(3*(tri_idx-1) + 2) = dead_bsearch_int(vert_gid, nverts_local, 2 + ix_r)
                tri_conn(3*(tri_idx-1) + 3) = dead_bsearch_int(vert_gid, nverts_local, 2 + i_elem)
                elem_gids(tri_idx) = elem_gid
                entity_order(tri_idx) = n
             else if (j_elem == nyg - 1) then
                ! North pole triangle: (i_elem, nyg-1), (ix_r, nyg-1), NP
                tri_idx = tri_idx + 1
                tri_conn(3*(tri_idx-1) + 1) = dead_bsearch_int(vert_gid, nverts_local, 2 + (nyg-2)*nxg + i_elem)
                tri_conn(3*(tri_idx-1) + 2) = dead_bsearch_int(vert_gid, nverts_local, 2 + (nyg-2)*nxg + ix_r)
                tri_conn(3*(tri_idx-1) + 3) = dead_bsearch_int(vert_gid, nverts_local, vgid_np)
                elem_gids(tri_idx) = elem_gid
                entity_order(tri_idx) = n
             else
                ! Interior quad: BL, BR, TR, TL (counter-clockwise)
                quad_idx = quad_idx + 1
                vgid_bl = 2 + (j_elem - 1) * nxg + i_elem
                vgid_br = 2 + (j_elem - 1) * nxg + ix_r
                vgid_tr = 2 + j_elem * nxg + ix_r
                vgid_tl = 2 + j_elem * nxg + i_elem
                quad_conn(4*(quad_idx-1) + 1) = dead_bsearch_int(vert_gid, nverts_local, vgid_bl)
                quad_conn(4*(quad_idx-1) + 2) = dead_bsearch_int(vert_gid, nverts_local, vgid_br)
                quad_conn(4*(quad_idx-1) + 3) = dead_bsearch_int(vert_gid, nverts_local, vgid_tr)
                quad_conn(4*(quad_idx-1) + 4) = dead_bsearch_int(vert_gid, nverts_local, vgid_tl)
                elem_gids(ntri_elems + quad_idx) = elem_gid
                entity_order(ntri_elems + quad_idx) = n
             endif
          enddo

          ! Create triangular elements (polar caps) - MBTRI = 2
          if (ntri_elems > 0) then
             elem_type = 2  ! MBTRI
             ierr = iMOAB_CreateElements(mbdomain, ntri_elems, elem_type, 3, tri_conn, block_ID=1)
             if (ierr .ne. 0) &
                call shr_sys_abort('Error: fail to create RLL triangle elements')
          endif

          ! Create quad elements (interior) - MBQUAD = 3
          if (nquad_elems > 0) then
             elem_type = 3  ! MBQUAD
             ierr = iMOAB_CreateElements(mbdomain, nquad_elems, elem_type, 4, quad_conn, block_ID=2)
             if (ierr .ne. 0) &
                call shr_sys_abort('Error: fail to create RLL quad elements')
          endif

          deallocate(tri_conn, quad_conn)

          ! Set GLOBAL_ID tag on vertices for parallel resolution
          tagname='GLOBAL_ID'//C_NULL_CHAR
          ierr = iMOAB_DefineTagStorage(mbdomain, tagname, 0, 1, tagindex)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to define GLOBAL_ID tag for vertices')

          ! Set GLOBAL_ID on vertices (etype=0)
          ierr = iMOAB_SetIntTagStorage(mbdomain, tagname, nverts_local, 0, vert_gid)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to set GLOBAL_ID tag for vertices')

          ! Set GLOBAL_ID on elements (etype=1)
          ! elem_gids is ordered: tris first, then quads (matches MOAB entity order)
          ierr = iMOAB_SetIntTagStorage(mbdomain, tagname, nelems, 1, elem_gids)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to set GLOBAL_ID tag for elements')

          ! Resolve shared entities across MPI ranks using global vertex IDs
          ierr = iMOAB_ResolveSharedEntities(mbdomain, nverts_local, vert_gid)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to resolve shared entities')

          deallocate(vert_coords, vert_gid, elem_gids)
       end block
    else

         ! Create MOAB mesh for dead domain
      allocate(moab_vert_coords(lsize*3))
      do n = 1,lsize
         lonv = gbuf(n,dead_grid_lon) * shr_const_pi / 180.0_R8
         latv = gbuf(n,dead_grid_lat) * shr_const_pi / 180.0_R8
         moab_vert_coords(3*n-2)=COS(latv)*COS(lonv)
         moab_vert_coords(3*n-1)=COS(latv)*SIN(lonv)
         moab_vert_coords(3*n  )=SIN(latv)
      enddo

      ! create the vertices with coordinates from MCT domain
      ierr = iMOAB_CreateVertices(mbdomain, lsize*3, 3, moab_vert_coords)
      if (ierr .ne. 0)  &
         call shr_sys_abort('Error: fail to create MOAB vertices in X-' // trim(model) // ' model')
      deallocate(moab_vert_coords)

      tagname='GLOBAL_ID'//C_NULL_CHAR
      ierr = iMOAB_DefineTagStorage(mbdomain, tagname, &
                                    0, & ! dense, integer
                                    1, & ! number of components
                                    tagindex )
      if (ierr .ne. 0)  &
         call shr_sys_abort('Error: fail to retrieve GLOBAL_ID tag ')

      allocate(idata(lsize))
      ! get list of global IDs for Dofs
      call mct_gsMap_orderedPoints(gsMap, my_task, idata)

      ierr = iMOAB_SetIntTagStorage ( mbdomain, tagname, lsize, &
                                       etype, & ! vertex type
                                       idata)
      if (ierr .ne. 0)  &
         call shr_sys_abort('Error: fail to set GLOBAL_ID tag ')

      deallocate(idata)

    end if

   !  if (fullmesh) then
   !    ierr = iMOAB_ResolveSharedEntities( mbdomain, lsize, idata );
   !    if (ierr .ne. 0)  &
   !       call shr_sys_abort('Error: fail to resolve shared entities')
   !  end if

    ierr = iMOAB_UpdateMeshInfo( mbdomain )
    if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to update mesh info ')

    ierr = iMOAB_DefineTagStorage( mbdomain, "lat:lon:area:aream:frac:mask"//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create tag: area:aream:frac:mask' )

    if (lsize .gt. 0) then

         call define_reset_fields_moab(mbdomain, model, x2d, d2x, mpicom, logunit )

         allocate(data(lsize))

         ! For fullmesh, reorder data to match MOAB entity iteration order
         ! (all MBTRI first, then all MBQUAD). entity_order maps MOAB
         ! element position -> local gbuf index.
         if (fullmesh) then
            do n = 1, lsize
               data(n) = gbuf(entity_order(n), dead_grid_lat)
            enddo
         else
            data(:) = gbuf(:,dead_grid_lat)
         endif
         tagname='lat'//C_NULL_CHAR
         ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                          etype, &
                                          data)
         if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to set lat tag ')

         if (fullmesh) then
            do n = 1, lsize
               data(n) = gbuf(entity_order(n), dead_grid_lon)
            enddo
         else
            data(:) = gbuf(:,dead_grid_lon)
         endif
         tagname='lon'//C_NULL_CHAR
         ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                          etype, &
                                          data)
         if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to set lon tag ')

         if (fullmesh) then
            do n = 1, lsize
               data(n) = gbuf(entity_order(n), dead_grid_area)
            enddo
         else
            data(:) = gbuf(:,dead_grid_area)
         endif
         tagname='area'//C_NULL_CHAR
         ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                          etype, &
                                          data)
         if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to get area tag ')

         ! set the same data for aream (model area) as area
         tagname='aream'//C_NULL_CHAR
         ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                          etype, &
                                          data)
         if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to set aream tag ')

         if (fullmesh) then
            do n = 1, lsize
               data(n) = gbuf(entity_order(n), dead_grid_mask)
            enddo
         else
            data(:) = gbuf(:,dead_grid_mask)
         endif
         tagname='mask'//C_NULL_CHAR
         ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                          etype, &
                                          data)
         if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to set mask tag ')

         if (fullmesh) then
            do n = 1, lsize
               data(n) = gbuf(entity_order(n), dead_grid_frac)
            enddo
         else
            data(:) = gbuf(:,dead_grid_frac)
         endif
         tagname='frac'//C_NULL_CHAR
         ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                          etype, &
                                          data)
         if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to set frac tag ')

         ! deallocate temporary arrays
         deallocate(data)
         if (allocated(entity_order)) deallocate(entity_order)

    end if


#ifdef MOABDEBUG
    !      debug test - write out the mesh file to disk
    wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mbdomain, trim(model) // '_unr_deadModels.h5m'//C_NULL_CHAR, trim(wopts))
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR in writing data mesh lnd ')
    endif
#endif

  end subroutine dead_init_moab
  !===============================================================================


  subroutine define_reset_fields_moab( mbdomain, model, x2d, d2x, mpicom, logunit )

    use iMOAB, only: iMOAB_DefineTagStorage, &
                     iMOAB_SetIntTagStorage, &
                     iMOAB_SetDoubleTagStorage, &
                     iMOAB_WriteMesh
      implicit none
    !-------------------------------------------------------------------
    !---arguments---
    integer(IN)    , intent(inout)  :: mbdomain
    character(len=*) , intent(in)    :: model
    type(mct_aVect)  , intent(in) :: x2d         ! driver -> dead
    type(mct_aVect)  , intent(in) :: d2x         ! dead   -> driver
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)    , intent(in)  :: logunit

    !---local variables---
    integer(IN)          :: ierr        ! error code
    integer(IN)          :: my_task     ! mpi task within communicator
    integer(IN)          :: lsize       ! temporary
    integer(IN)          :: n    ,j,i   ! indices
    integer(IN)          :: ier         ! error status
    character(CL)        :: tagname
    integer(IN)          :: tagindex, etype
    integer                 :: nflds_d2x, nflds_x2d
    character(len=2048)    :: flds_x2d
    character(len=2048)    :: flds_d2x
    real(R8), pointer    :: data_x2d(:,:), data_d2x(:,:)     ! temporary

#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
    character(*), parameter :: subName = "(define_reset_fields_moab) "
    !-------------------------------------------------------------------
    !
    ! Initialize MOAB dead domain
    !

    ! Allocate memory
    !
    lsize = mct_avect_lsize(x2d)

    flds_d2x = trim(mct_aVect_exportRList2c(d2x)) ! exportRListToChar
    nflds_d2x = mct_avect_nRattr(d2x)

    flds_x2d = trim(mct_aVect_exportRList2c(x2d)) ! exportRListToChar
    nflds_x2d = mct_avect_nRattr(x2d)

    !
    ! Initialize attribute vector with special value
    !
    call mpi_comm_rank(mpicom, my_task, ier)

    etype = 1  ! cells (all dead components use element-based tags)

    ! first define the tags for x2d and d2x fields
    ierr = iMOAB_DefineTagStorage( mbdomain, trim(flds_x2d)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create tag: ' // trim(flds_x2d) )

    ierr = iMOAB_DefineTagStorage( mbdomain, trim(flds_d2x)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create tag: ' // trim(flds_d2x) )

    if (lsize .gt. 0) then
       allocate(data_x2d(lsize,nflds_x2d), data_d2x(lsize,nflds_d2x))

       !  data(:) = -9999.0_R8  ! generic special value
       data_x2d(:,:) = 0.0_R8  ! generic special value
       data_d2x(:,:) = 0.0_R8  ! generic special value

       ierr = iMOAB_SetDoubleTagStorage( mbdomain, trim(flds_x2d)//C_NULL_CHAR, &
                                       lsize*nflds_x2d, & ! size
                                       etype, & ! set data on vertices
                                       data_x2d )
       if (ierr > 0 )  &
          call shr_sys_abort('Error: fail to reset tag: ' // trim(flds_x2d) )


       ierr = iMOAB_SetDoubleTagStorage( mbdomain, trim(flds_d2x)//C_NULL_CHAR, &
                                          lsize*nflds_d2x, & ! size
                                          etype, & ! set data on vertices
                                          data_d2x )
       if (ierr > 0 )  &
            call shr_sys_abort('Error: fail to reset tag: ' // trim(flds_d2x) )

       deallocate(data_x2d)
       deallocate(data_d2x)

     end if

#ifdef MOABDEBUG
       !      debug test
    outfile = trim(model) // '_deadModels.h5m'//C_NULL_CHAR
    wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
       !      write out the mesh file to disk
    ierr = iMOAB_WriteMesh(mbdomain, trim(outfile), trim(wopts))
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR in writing data mesh lnd ')
    endif
#endif

  end subroutine define_reset_fields_moab
  !===============================================================================

  !-----------------------------------------------------------------------
  ! Heapsort for integer arrays. Sorts arr(1:n) in-place in ascending order.
  ! O(n log n) worst case, no extra memory.
  !-----------------------------------------------------------------------
  subroutine dead_heapsort_int(arr, n)
    use shr_kind_mod, only : IN=>SHR_KIND_IN
    implicit none
    integer(IN), intent(in)    :: n
    integer(IN), intent(inout) :: arr(n)
    integer(IN) :: i, ir, j, l, tmp

    if (n < 2) return

    ! Build max-heap
    l = n/2 + 1
    ir = n
    do
       if (l > 1) then
          l = l - 1
          tmp = arr(l)
       else
          tmp = arr(ir)
          arr(ir) = arr(1)
          ir = ir - 1
          if (ir == 1) then
             arr(1) = tmp
             return
          endif
       endif
       i = l
       j = l + l
       do while (j <= ir)
          if (j < ir) then
             if (arr(j) < arr(j+1)) j = j + 1
          endif
          if (tmp < arr(j)) then
             arr(i) = arr(j)
             i = j
             j = j + j
          else
             j = ir + 1
          endif
       enddo
       arr(i) = tmp
    enddo
  end subroutine dead_heapsort_int

  !-----------------------------------------------------------------------
  ! Binary search for integer value in sorted array arr(1:n).
  ! Returns the 1-based index of val. Aborts if val is not found.
  !-----------------------------------------------------------------------
  pure function dead_bsearch_int(arr, n, val) result(idx)
    use shr_kind_mod, only : IN=>SHR_KIND_IN
    implicit none
    integer(IN), intent(in) :: n
    integer(IN), intent(in) :: arr(n)
    integer(IN), intent(in) :: val
    integer(IN) :: idx
    integer(IN) :: lo, hi, mid

    lo = 1
    hi = n
    do while (lo <= hi)
       mid = (lo + hi) / 2
       if (arr(mid) < val) then
          lo = mid + 1
       else if (arr(mid) > val) then
          hi = mid - 1
       else
          idx = mid
          return
       endif
    enddo
    ! Should never reach here if mesh construction is correct
    idx = -1
  end function dead_bsearch_int

#endif

end module dead_mct_mod
