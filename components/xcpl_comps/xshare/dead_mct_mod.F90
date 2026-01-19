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
#endif

!===============================================================================
contains
!===============================================================================

  subroutine dead_init_mct(model, Eclock, x2d, d2x, &
         flds_x2d, flds_d2x, &
         gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, nxg, nyg )

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
    logical                  :: flood=.false. ! rof flood flag

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
       nxg, nyg, decomp_type, nproc_x, seg_len, flood)

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
    integer(IN)      :: ent_type = 0    ! entity type: default = vertex
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

    lsize = mct_avect_lsize(x2d)
    nflds_d2x = mct_avect_nRattr(d2x)
    nflds_x2d = mct_avect_nRattr(x2d)

#ifdef HAVE_MOAB
    flds_d2x = trim(mct_aVect_exportRList2c(d2x)) ! exportRListToChar
    flds_x2d = trim(mct_aVect_exportRList2c(x2d)) ! exportRListToChar

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
         do n=1,lsize
           data_d2x(n,ki) = d2x%rAttr(ki,n)
         enddo
#endif

    case('glc')

       ki = mct_aVect_indexRA(d2x,"Sg_icemask",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8
#ifdef HAVE_MOAB
       data_d2x(:,ki) = 1.0_R8
#endif

       ki = mct_aVect_indexRA(d2x,"Sg_icemask_coupled_fluxes",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8
#ifdef HAVE_MOAB
       data_d2x(:,ki) = 1.0_R8
#endif

       ki = mct_aVect_indexRA(d2x,"Sg_ice_covered",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8
#ifdef HAVE_MOAB
       data_d2x(:,ki) = 1.0_R8
#endif
    end select

#ifdef HAVE_MOAB
    ! Import data into MOAB structures
   !  call define_import_fields_moab(model, mpicom, gbuf, gsMap, logunit, ggrid, &
   !       flds_x2d, flds_d2x, data_x2d, data_d2x, nxg, nyg)

   !  print *, ' DEAD RUN MCT: setting MOAB tags for model ', trim(model), ' nflds_d2x=', nflds_d2x, ' lsize=', lsize
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, trim(flds_d2x) // C_NULL_CHAR, lsize*nflds_d2x, &
                                     ent_type, & ! NOTE: setting data on vertices
                                     data_d2x)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set flds_d2x fields tag ', ierr)

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
  subroutine dead_init_moab( mbdomain, model, gsMap, gbuf, flds_x2d, flds_d2x, mpicom, compid, logunit, nxg, nyg )

    use iMOAB, only: iMOAB_CreateVertices, &
                     iMOAB_CreateElements, &
                     iMOAB_DefineTagStorage, &
                     iMOAB_SetIntTagStorage, &
                     iMOAB_SetDoubleTagStorage, &
                     iMOAB_ResolveSharedEntities, &
                     iMOAB_UpdateMeshInfo, &
                     iMOAB_GetMeshInfo, &
                     iMOAB_WriteMesh
      implicit none
    !-------------------------------------------------------------------
    !---arguments---
    integer(IN)    , intent(in)   :: mbdomain
    character(len=*) , intent(in) :: model
    type(mct_gsMap), intent(in)   :: gsMap
    real(R8)       , intent(in)   :: gbuf(:,:)
    character(len=*) , intent(in)    :: flds_x2d
    character(len=*) , intent(in)    :: flds_d2x
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
    logical          :: fullmesh =.false.
! #ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
! #endif
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

    etype = 0  ! vertices
   !  if (trim(model) .ne. 'atm' .and. trim(model) .ne. 'lnd' .and. trim(model) .ne. 'rof') then
   !    fullmesh = .true.
   !    etype = 1  ! elements
   !  end if

    if (fullmesh) then
       ! Create RLL dual mesh: centroids (lat,lon) become element centers; vertices at corners
       ! For nxg x nyg centroid grid, we need (nxg+1) x (nyg+1) corner vertices
       ! Each centroid becomes a quad element using its 4 neighboring corners
       block
          integer(IN)          :: nverts_total   ! total vertices in RLL mesh
          integer(IN)          :: nelems         ! number of elements
          integer(IN), allocatable :: elem_conn(:,:)  ! element connectivity
          integer(IN), allocatable :: vert_glid(:)    ! global IDs for vertices
          integer(IN)          :: ie, iv, ix, iy, v0, v1, v2, v3
          integer(IN)          :: elem_type
          integer(IN), allocatable :: elem_conn_flat(:)
          real(R8), allocatable :: vert_coords(:)     ! vertex coords in 3D
          real(R8)             :: lat_c, lon_c, lat_v, lon_v
          real(R8)             :: dlat, dlon          ! grid spacing
          integer(IN)          :: i_global, j_global  ! global grid indices for centroid
          logical, allocatable :: vert_in_local(:)    ! which vertices belong locally

          ! Grid spacing (approximate; valid for regular grid)
          dlat = 180.0_R8 / real(nyg, R8)
          dlon = 360.0_R8 / real(nxg, R8)

          ! Vertices in full RLL mesh (corners only exist if centroid exists locally or is adjacent)
          nverts_total = (nxg + 1) * (nyg + 1)
          nelems = nxg * nyg

          allocate(vert_glid(nverts_total))
          allocate(vert_in_local(nverts_total))
          allocate(vert_coords(nverts_total * 3))
          allocate(elem_conn(4, nelems))
          allocate(elem_conn_flat(nelems * 4))

          ! Initialize: mark all vertices as non-local initially
          vert_in_local = .false.
          vert_coords = 0.0_R8

          ! Loop through local centroid points and mark adjacent vertices
          do n = 1, lsize
             lat_c = gbuf(n, dead_grid_lat)
             lon_c = gbuf(n, dead_grid_lon)
             ! Find global grid position (i_global, j_global) of this centroid in [0, nxg-1] x [0, nyg-1]
             ! Use MCT global index to infer position
             ! For simplicity, assume sequential ordering: centroid (i,j) has global index i + j*nxg
             i_global = mod(n-1, nxg)
             j_global = (n-1) / nxg

             ! Mark 4 corner vertices of this centroid's cell as local
             ! Vertices are in [0, nxg] x [0, nyg], so corners of cell (i,j) are:
             !   (i, j), (i+1, j), (i+1, j+1), (i, j+1)
             do iy = j_global, j_global + 1
                do ix = i_global, i_global + 1
                   if (ix >= 0 .and. ix <= nxg .and. iy >= 0 .and. iy <= nyg) then
                      iv = iy * (nxg + 1) + ix + 1  ! 1-indexed vertex ID
                      vert_in_local(iv) = .true.
                      ! Compute corner lat/lon (offset by half-grid spacing from centroid)
                      lat_v = lat_c + (real(iy, R8) - real(j_global, R8) - 0.5_R8) * dlat
                      lon_v = lon_c + (real(ix, R8) - real(i_global, R8) - 0.5_R8) * dlon
                      ! Clamp lat to [-90, 90]
                      lat_v = max(-90.0_R8, min(90.0_R8, lat_v))
                      ! Convert to 3D spherical coords
                      vert_coords(3*iv-2) = cos(lat_v * shr_const_pi / 180.0_R8) * cos(lon_v * shr_const_pi / 180.0_R8)
                      vert_coords(3*iv-1) = cos(lat_v * shr_const_pi / 180.0_R8) * sin(lon_v * shr_const_pi / 180.0_R8)
                      vert_coords(3*iv  ) = sin(lat_v * shr_const_pi / 180.0_R8)
                      ! Global vertex ID = global grid position
                      vert_glid(iv) = iy * (nxg + 1) + ix + 1
                   endif
                enddo
             enddo
          enddo

          ! Create vertices for local mesh (only the ones marked as local)
          allocate(moab_vert_coords(count(vert_in_local)*3))
          iv = 0
          do i = 1, nverts_total
             if (vert_in_local(i)) then
                iv = iv + 1
                moab_vert_coords(3*iv-2:3*iv) = vert_coords(3*i-2:3*i)
             endif
          enddo
          ierr = iMOAB_CreateVertices(mbdomain, iv*3, 3, moab_vert_coords)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to create corner vertices for RLL mesh')
          deallocate(moab_vert_coords)

          ! Create quad elements: each centroid becomes a quad using its 4 corner vertices
          ie = 0
          do n = 1, lsize
             i_global = mod(n-1, nxg)
             j_global = (n-1) / nxg
             ! 4 corners of element (i,j): at vertex positions (i,j), (i+1,j), (i+1,j+1), (i,j+1)
             v0 = j_global * (nxg + 1) + i_global + 1      ! (i, j)
             v1 = j_global * (nxg + 1) + i_global + 1 + 1  ! (i+1, j)
             v2 = (j_global+1) * (nxg + 1) + i_global + 1 + 1  ! (i+1, j+1)
             v3 = (j_global+1) * (nxg + 1) + i_global + 1  ! (i, j+1)

             ie = ie + 1
             elem_conn(1, ie) = v0 - 1  ! Convert to 0-indexed
             elem_conn(2, ie) = v1 - 1
             elem_conn(3, ie) = v2 - 1
             elem_conn(4, ie) = v3 - 1
          enddo

          ! Flatten connectivity
          elem_conn_flat = reshape(elem_conn, [nelems * 4])

          ! Create elements
          elem_type = 3  ! quad
          ierr = iMOAB_CreateElements(mbdomain, lsize, elem_type, 4, elem_conn_flat, block_ID=1)
          if (ierr .ne. 0) &
             call shr_sys_abort('Error: fail to create RLL quad elements')

          deallocate(elem_conn)
          deallocate(elem_conn_flat)
          deallocate(vert_coords)
          deallocate(vert_glid)
          deallocate(vert_in_local)
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
         call shr_sys_abort('Error: fail to create MOAB vertices in data lnd model')
      deallocate(moab_vert_coords)

    end if

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
                                     0, & ! vertex type
                                     idata)
    if (ierr .ne. 0)  &
       call shr_sys_abort('Error: fail to set GLOBAL_ID tag ')

    if (fullmesh) then
      ierr = iMOAB_ResolveSharedEntities( mbdomain, lsize, idata );
      if (ierr .ne. 0)  &
         call shr_sys_abort('Error: fail to resolve shared entities')
    end if

    ierr = iMOAB_UpdateMeshInfo( mbdomain )
    if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to update mesh info ')

    deallocate(idata)

    allocate(data(lsize))
    ierr = iMOAB_DefineTagStorage( mbdomain, "lat:lon:area:aream:frac:mask"//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to create tag: area:aream:frac:mask' )


    data(:) = gbuf(:,dead_grid_lat)
    tagname='lat'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                    0, & ! set data on vertices
                                    data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set lat tag ')

    data(:) = gbuf(:,dead_grid_lon)
    tagname='lon'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                    0, & ! set data on vertices
                                    data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set lon tag ')

    data(:) = gbuf(:,dead_grid_area)
    tagname='area'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                     etype, & ! set data on elements
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to get area tag ')

    ! set the same data for aream (model area) as area
    ! data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'aream'),:)
    tagname='aream'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                     etype, & ! set data on elements
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set aream tag ')

    data(:) = gbuf(:,dead_grid_mask)
    tagname='mask'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                     etype, & ! set data on elements
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set mask tag ')

    data(:) = gbuf(:,dead_grid_frac)
    tagname='frac'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
                                     etype, & ! set data on elements
                                     data)
    if (ierr > 0 )  &
       call shr_sys_abort('Error: fail to set frac tag ')

    ! deallocate temporary arrays
    deallocate(data)

    call define_reset_fields_moab(mbdomain, model, gbuf, flds_x2d, flds_d2x, mpicom, logunit, nxg, nyg )

  !      debug test
   !  outfile = trim(model) // '_deadModels.h5m'//C_NULL_CHAR
    wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
   !  wopts   = ''//C_NULL_CHAR !
       !      write out the mesh file to disk
    ierr = iMOAB_WriteMesh(mbdomain, trim(model) // '_unr_deadModels.h5m'//C_NULL_CHAR, trim(wopts))
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR in writing data mesh lnd ')
    endif

  end subroutine dead_init_moab
  !===============================================================================


  subroutine define_reset_fields_moab( mbdomain, model, gbuf, flds_x2d, flds_d2x, mpicom, logunit, nxg, nyg )

    use iMOAB, only: iMOAB_DefineTagStorage, &
                     iMOAB_SetIntTagStorage, &
                     iMOAB_SetDoubleTagStorage, &
                     iMOAB_WriteMesh
      implicit none
    !-------------------------------------------------------------------
    !---arguments---
    integer(IN)    , intent(in)  :: mbdomain
    character(len=*) , intent(in)    :: model
    real(R8)       , intent(in)  :: gbuf(:,:)
    character(len=*) , intent(in)    :: flds_x2d
    character(len=*) , intent(in)    :: flds_d2x
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)    , intent(in)  :: logunit
    integer(IN)    , intent(in)  :: nxg          ! global grid size in x-direction
    integer(IN)    , intent(in)  :: nyg          ! global grid size in y-direction

    !---local variables---
    integer(IN)          :: ierr        ! error code
    integer(IN)          :: my_task     ! mpi task within communicator
    integer(IN)          :: lsize       ! temporary
    integer(IN)          :: n    ,j,i   ! indices
    integer(IN)          :: ier         ! error status
    character(CL)        :: tagname
    integer(IN)          :: tagindex, etype
    real(R8), pointer    :: data(:)     ! temporary

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
    lsize = size(gbuf,dim=1)
    allocate(data(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mpi_comm_rank(mpicom, my_task, ier)
   !  call mct_aVect_zero(domain%data)
    data(:) = -9999.0_R8  ! generic special value

    data(:) = 0.0_R8  ! generic special value

    etype = 0  ! vertices
   !  if (trim(model) .ne. 'atm' .and. trim(model) .ne. 'lnd' .and. trim(model) .ne. 'rof') then
   !    fullmesh = .true.
   !    etype = 1  ! cells
   !  end if

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
       call shr_sys_abort('Error: fail to create tag: ' // trim(flds_x2d) )

   !  ! Now let us set the default values for these fields to zero
   !  data(:) = gbuf(:,dead_grid_lat)
   !  tagname='lat'//C_NULL_CHAR
   !  ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
   !                                  etype, & ! set data on vertices
   !                                  data)
   !  if (ierr > 0 )  &
   !     call shr_sys_abort('Error: fail to set lat tag ')

   !  data(:) = gbuf(:,dead_grid_lon)
   !  tagname='lon'//C_NULL_CHAR
   !  ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
   !                                  0, & ! set data on vertices
   !                                  data)
   !  if (ierr > 0 )  &
   !     call shr_sys_abort('Error: fail to set lon tag ')

   !  data(:) = gbuf(:,dead_grid_area)
   !  tagname='area'//C_NULL_CHAR
   !  ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
   !                                   etype, & ! set data on elements
   !                                   data)
   !  if (ierr > 0 )  &
   !     call shr_sys_abort('Error: fail to get area tag ')

   !  ! set the same data for aream (model area) as area
   !  ! data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'aream'),:)
   !  tagname='aream'//C_NULL_CHAR
   !  ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
   !                                   etype, & ! set data on elements
   !                                   data)
   !  if (ierr > 0 )  &
   !     call shr_sys_abort('Error: fail to set aream tag ')

   !  data(:) = gbuf(:,dead_grid_mask)
   !  tagname='mask'//C_NULL_CHAR
   !  ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
   !                                   etype, & ! set data on elements
   !                                   data)
   !  if (ierr > 0 )  &
   !     call shr_sys_abort('Error: fail to set mask tag ')

   !  data(:) = gbuf(:,dead_grid_frac)
   !  tagname='frac'//C_NULL_CHAR
   !  ierr = iMOAB_SetDoubleTagStorage ( mbdomain, tagname, lsize, &
   !                                   etype, & ! set data on elements
   !                                   data)
   !  if (ierr > 0 )  &
   !     call shr_sys_abort('Error: fail to set frac tag ')

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

    deallocate(data)

  end subroutine define_reset_fields_moab
  !===============================================================================

#endif

end module dead_mct_mod
