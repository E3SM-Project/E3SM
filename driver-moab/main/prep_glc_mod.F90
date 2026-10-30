module prep_glc_mod

#include "shr_assert.h"
  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_kind_mod    , only: CXX => SHR_KIND_CXX
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_glc, num_inst_lnd, num_inst_frc, &
                              num_inst_ocn
  use seq_comm_mct    , only: CPLID, GLCID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs
  use seq_comm_mct    , only: mblxid, mbgxid, mbintxlg, mbintxgl ! iMOAB app ids: lnd and glc on coupler, l2g and g2l map holders
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: glc, lnd, ocn
  use glc_elevclass_mod, only : glc_get_num_elevation_classes, glc_elevclass_as_string

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_glc_init
  public :: prep_glc_mrg_ocn
  public :: prep_glc_mrg_lnd

  public :: prep_glc_accum_ocn
  public :: prep_glc_accum_lnd
  public :: prep_glc_accum_avg

  public :: prep_glc_calc_o2x_gx
  public :: prep_glc_calc_l2x_gx

  public :: prep_glc_zero_fields

  public :: prep_glc_get_l2gacc_lm
  public :: prep_glc_get_l2gacc_lm_cnt
  public :: prep_glc_get_l2gacc_lm_cnt_avg
  public :: prep_glc_get_sharedFieldsLndGlc

  public :: prep_glc_get_l2x_gx
  public :: prep_glc_get_l2gacc_lx
  public :: prep_glc_get_l2gacc_lx_one_instance
  public :: prep_glc_get_l2gacc_lx_cnt
  public :: prep_glc_get_l2gacc_lx_cnt_avg

  public :: prep_glc_get_o2x_gx
  public :: prep_glc_get_x2gacc_gx
  public :: prep_glc_get_x2gacc_gx_cnt

  public :: prep_glc_get_mapper_Sl2g
  public :: prep_glc_get_mapper_Fl2g

  public :: prep_glc_get_mapper_So2g_shelf
  public :: prep_glc_get_mapper_Fo2g_shelf

  public :: prep_glc_calculate_subshelf_boundary_fluxes

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_do_renormalize_smb
  private :: prep_glc_renormalize_smb

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sl2g
  type(seq_map), pointer :: mapper_Fl2g
  type(seq_map), pointer :: mapper_So2g_shelf
  type(seq_map), pointer :: mapper_Fo2g_shelf
  type(seq_map), pointer :: mapper_So2g_tf
  type(seq_map), pointer :: mapper_Fg2l

  ! attribute vectors
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: o2x_gx(:) ! Ocn export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  
  type(mct_aVect), pointer :: x2gacc_gx(:) ! Glc export, glc grid, cpl pes - allocated in driver
  integer        , target :: x2gacc_gx_cnt ! x2gacc_gx: number of time samples accumulated 

  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated
  integer        , target :: l2gacc_lx_cnt_avg ! l2gacc_lx: number of time samples averaged

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator

  ! Whether to renormalize the SMB for conservation.
  ! Should be set to true for 2-way coupled runs with evolving ice sheets.
  ! Does not need to be true for 1-way coupling.
  logical :: smb_renormalize
  logical :: glc_present ! .true. => glc is present

  ! Name of flux field giving surface mass balance
  character(len=*), parameter :: qice_fieldname = 'Flgl_qice'

  ! Names of some other fields
  character(len=*), parameter :: Sg_frac_field = 'Sg_ice_covered'
  character(len=*), parameter :: Sg_topo_field = 'Sg_topo'
  character(len=*), parameter :: Sg_icemask_field = 'Sg_icemask'

  type(mct_aVect), pointer :: o2gacc_ox(:) ! Ocn export, lnd grid, cpl pes - allocated in driver
  integer        , target :: o2gacc_ox_cnt ! number of time samples accumulated

  real(r8), allocatable ::  oceanTemperature(:)
  real(r8), allocatable ::  oceanSalinity(:)
  real(r8), allocatable ::  oceanHeatTransferVelocity(:)
  real(r8), allocatable ::  oceanSaltTransferVelocity(:)
  real(r8), allocatable ::  interfacePressure(:)
  real(r8), allocatable ::  iceTemperature(:)
  real(r8), allocatable ::  iceTemperatureDistance(:)
  integer, allocatable ::  iceFloatingMask(:)
  real(r8), allocatable ::  outInterfaceSalinity(:)
  real(r8), allocatable ::  outInterfaceTemperature(:)
  real(r8), allocatable ::  outFreshwaterFlux(:)
  real(r8), allocatable ::  outOceanHeatFlux(:)
  real(r8), allocatable ::  outIceHeatFlux(:)

  ! moab support: accumulation of the per-elevation-class lnd fields on the coupler
  ! land mesh, and bookkeeping for the lnd->glc downscaling done on moab tags
  character(CXX)        :: sharedFieldsLndGlc      ! = seq_flds_l2x_fields_to_glc, the tag list
  real(r8), allocatable, target :: l2gacc_lm(:,:)  ! accumulated lnd fields, (lsize_lm, nflds_lg)
  real(r8), allocatable :: l2x_lm2(:,:)            ! scratch for reading the instantaneous lnd tags
  integer , target      :: l2gacc_lm_cnt           ! l2gacc_lm: number of time samples accumulated
  integer , target      :: l2gacc_lm_cnt_avg       ! l2gacc_lm: number of samples in last average
  integer               :: lsize_lm = 0            ! number of local cells, land coupler mesh
  integer               :: lsize_gm = 0            ! number of local cells, glc coupler mesh
  integer               :: nflds_lg = 0            ! number of fields in sharedFieldsLndGlc

  !================================================================================================

contains

  !================================================================================================

  subroutine prep_glc_init(infodata, lnd_c2_glc, ocn_c2_glctf, ocn_c2_glcshelf)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and mapping variables
    !
    use iMOAB, only : iMOAB_RegisterApplication, iMOAB_MigrateMapMesh, &
         iMOAB_ComputeCommGraph, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo
    use iso_c_binding, only : C_NULL_CHAR
    use shr_string_mod, only : shr_string_listGetNum
    use map_glc2lnd_mod, only : map_glc2lnd_ec_init
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: lnd_c2_glc ! .true.  => lnd to glc coupling on
    logical                  , intent(in)    :: ocn_c2_glctf ! .true.  => ocn to glc thermal forcing coupling on
    logical                  , intent(in)    :: ocn_c2_glcshelf ! .true.  => ocn to glc shelf coupling on
    !
    ! Local Variables
    integer                          :: eli, egi, eoi
    integer                          :: lsize_l
    integer                          :: lsize_g
    integer                         ::  lsize_o
    logical                          :: samegrid_lg   ! samegrid land and glc
    logical                          :: samegrid_go   ! .true. => samegrid ocean and glc
    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                          :: do_hist_l2x1yrg ! .true. => create aux files: l2x 1yr glc forcings
    character(CL)                    :: lnd_gnam      ! lnd grid
    character(CL)                    :: glc_gnam      ! glc grid
    character(CL)                    :: ocn_gnam      ! ocn grid

    type(mct_avect), pointer         :: l2x_lx
    type(mct_avect), pointer         :: x2g_gx
    type(mct_avect), pointer         :: o2x_ox

    ! moab locals
    integer                          :: ierr, idintx, type_grid, arearead
    integer                          :: tagtype, numco, tagindex
    integer                          :: mpigrp_CPLID  ! coupler pes group
    integer                          :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
    character(CL)                    :: appname
    character(CL)                    :: wgtIdSl2g, wgtIdFl2g, wgtIdFg2l
    character(CXX)                   :: tagname

    character(*), parameter          :: subname = '(prep_glc_init)'
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         glc_present=glc_present       , &
         histaux_l2x1yrg=do_hist_l2x1yrg, &
         lnd_gnam=lnd_gnam             , &
         glc_gnam=glc_gnam             , &
         ocn_gnam=ocn_gnam)

    allocate(mapper_Sl2g)
    allocate(mapper_Fl2g)
    allocate(mapper_So2g_shelf)
    allocate(mapper_So2g_tf)
    allocate(mapper_Fo2g_shelf)
    allocate(mapper_Fg2l)

    smb_renormalize = prep_glc_do_renormalize_smb(infodata)

    if ((glc_present .and. lnd_c2_glc) .or. do_hist_l2x1yrg) then

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       allocate(l2gacc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_init(l2gacc_lx(eli), rList=seq_flds_l2x_fields_to_glc, lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       end do
       l2gacc_lx_cnt = 0
       l2gacc_lx_cnt_avg = 0

       ! moab accumulator: the per-elevation-class lnd fields are accumulated in a
       ! plain array read from the coupler land mesh tags (prep_rof pattern)
       if (mblxid >= 0) then
          sharedFieldsLndGlc = trim(seq_flds_l2x_fields_to_glc)
          nflds_lg = shr_string_listGetNum(sharedFieldsLndGlc)
          ierr = iMOAB_GetMeshInfo ( mblxid, nvert, nvise, nbl, nsurf, nvisBC )
          if (ierr .ne. 0) then
             call shr_sys_abort(subname//' ERROR getting land coupler mesh info')
          endif
          lsize_lm = nvise(1)
          allocate(l2gacc_lm(lsize_lm, nflds_lg))
          allocate(l2x_lm2(lsize_lm, nflds_lg))
          l2gacc_lm = 0._r8
          l2x_lm2 = 0._r8
          l2gacc_lm_cnt = 0
          l2gacc_lm_cnt_avg = 0
       endif
    end if

    if (glc_present .and. lnd_c2_glc) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)

       allocate(l2x_gx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_init(l2x_gx(eli), rList=seq_flds_x2g_fields, lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))
       enddo

       if (lnd_c2_glc) then

          samegrid_lg = .true.
          if (trim(lnd_gnam) /= trim(glc_gnam)) samegrid_lg = .false.

          ! the mct sparse-matrix init (seq_map_init_rcfile) is not usable in
          ! driver-moab (coupler-side gsmaps are not populated); the mappers get
          ! their real (moab) context below
          call seq_map_mapinit(mapper_Sl2g, mpicom_CPLID)
          call seq_map_mapinit(mapper_Fl2g, mpicom_CPLID)
          ! We need to initialize our own Fg2l mapper because in some cases (particularly
          ! TG compsets - dlnd forcing CISM) the system doesn't otherwise create a Fg2l
          ! mapper.
          call seq_map_mapinit(mapper_Fg2l, mpicom_CPLID)

          ! now give the mappers moab context and load the mapping weights with iMOAB
          if ((mblxid >= 0) .and. (mbgxid >= 0)) then

             if (samegrid_lg) then
                call shr_sys_abort(subname// &
                     ' ERROR: moab lnd->glc coupling requires distinct lnd and glc grids with map files')
             endif

             call seq_comm_getData(CPLID, mpigrp=mpigrp_CPLID)
             type_grid = 3 ! FV-FV for both lnd and glc coupler meshes

             if (iamroot_CPLID) then
                write(logunit,*) ' '
                write(logunit,F00) 'Initializing MOAB mapper_Sl2g and mapper_Fl2g'
             end if
             appname = "LND_GLC_COU"//C_NULL_CHAR
             ! unique external id for the moab app holding the lnd->glc read map
             idintx = 100*lnd(1)%cplcompid + glc(1)%cplcompid
             ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxlg)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error in registering lnd glc map app'
                call shr_sys_abort(subname//' ERROR in registering lnd glc map app')
             endif

             ! scalar (bilinear) map, used for all the lnd->glc downscaling maps
             mapper_Sl2g%src_mbid = mblxid
             mapper_Sl2g%tgt_mbid = mbgxid
             mapper_Sl2g%intx_mbid = mbintxlg
             mapper_Sl2g%src_context = lnd(1)%cplcompid
             mapper_Sl2g%intx_context = idintx
             wgtIdSl2g = 'scalar_l2g'
             mapper_Sl2g%weight_identifier = wgtIdSl2g
             mapper_Sl2g%mbname = 'mapper_Sl2g'
             arearead = 0 ! no need for areas
             call moab_map_init_rcfile( mapper_Sl2g, type_grid, &
                  'seq_maps.rc', 'lnd2glc_smapname:', 'lnd2glc_smaptype:', samegrid_lg, &
                  arearead, wgtIdSl2g, 'mapper_Sl2g MOAB initialization', esmf_map_flag)

             ! flux (conservative) map; arearead=2 loads the map file area_b into the glc
             ! mesh aream, reproducing the mct driver's seq_map_readdata of area_b -> aream
             mapper_Fl2g%src_mbid = mblxid
             mapper_Fl2g%tgt_mbid = mbgxid
             mapper_Fl2g%intx_mbid = mbintxlg
             mapper_Fl2g%src_context = lnd(1)%cplcompid
             mapper_Fl2g%intx_context = idintx
             wgtIdFl2g = 'flux_l2g'
             mapper_Fl2g%weight_identifier = wgtIdFl2g
             mapper_Fl2g%mbname = 'mapper_Fl2g'
             arearead = 2 ! area_b for glc aream
             call moab_map_init_rcfile( mapper_Fl2g, type_grid, &
                  'seq_maps.rc', 'lnd2glc_fmapname:', 'lnd2glc_fmaptype:', samegrid_lg, &
                  arearead, wgtIdFl2g, 'mapper_Fl2g MOAB initialization', esmf_map_flag)

             ! one mesh migration and one comm graph cover both maps (same coverage)
             ierr = iMOAB_MigrateMapMesh (mblxid, mbintxlg, mpicom_CPLID, mpigrp_CPLID, &
                  mpigrp_CPLID, type_grid, lnd(1)%cplcompid, idintx)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error in migrating lnd mesh for map lnd 2 glc'
                call shr_sys_abort(subname//' ERROR in migrating lnd mesh for map lnd 2 glc')
             endif
             ierr = iMOAB_ComputeCommGraph( mblxid, mbintxlg, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, &
                  type_grid, type_grid, lnd(1)%cplcompid, idintx)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing comm graph for second hop, LND-GLC'
                call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, LND-GLC')
             endif

             ! define the per-elevation-class projection target tags on the glc mesh
             tagtype = 1  ! dense, double
             numco = 1
             tagname = trim(seq_flds_l2x_fields_to_glc)//C_NULL_CHAR
             ierr = iMOAB_DefineTagStorage(mbgxid, tagname, tagtype, numco, tagindex )
             if (ierr .ne. 0) then
                call shr_sys_abort(subname//' ERROR defining lnd per-EC tags on the glc mesh')
             endif

             ! glc coupler mesh size, used by the downscaling
             ierr = iMOAB_GetMeshInfo ( mbgxid, nvert, nvise, nbl, nsurf, nvisBC )
             if (ierr .ne. 0) then
                call shr_sys_abort(subname//' ERROR getting glc coupler mesh info')
             endif
             lsize_gm = nvise(1)

             ! moab context for the glc->lnd conservative map (used here by the smb
             ! renormalization). If prep_lnd has set it up already (glc_c2_lnd), just
             ! point our mapper at the same map app and weights; otherwise (e.g. TG
             ! compsets) register and load it ourselves.
             wgtIdFg2l = 'flux_g2l'
             if (mbintxgl < 0) then
                if (iamroot_CPLID) then
                   write(logunit,*) ' '
                   write(logunit,F00) 'Initializing MOAB mapper_Fg2l (glc->lnd map app)'
                end if
                appname = "GLC_LND_COU"//C_NULL_CHAR
                idintx = 100*glc(1)%cplcompid + lnd(1)%cplcompid
                ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxgl)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error in registering glc lnd map app'
                   call shr_sys_abort(subname//' ERROR in registering glc lnd map app')
                endif
                mapper_Fg2l%src_mbid = mbgxid
                mapper_Fg2l%tgt_mbid = mblxid
                mapper_Fg2l%intx_mbid = mbintxgl
                mapper_Fg2l%src_context = glc(1)%cplcompid
                mapper_Fg2l%intx_context = idintx
                mapper_Fg2l%weight_identifier = wgtIdFg2l
                mapper_Fg2l%mbname = 'mapper_Fg2l'
                arearead = 0
                call moab_map_init_rcfile( mapper_Fg2l, type_grid, &
                     'seq_maps.rc', 'glc2lnd_fmapname:', 'glc2lnd_fmaptype:', samegrid_lg, &
                     arearead, wgtIdFg2l, 'mapper_Fg2l (prep_glc) MOAB initialization', esmf_map_flag)
                ierr = iMOAB_MigrateMapMesh (mbgxid, mbintxgl, mpicom_CPLID, mpigrp_CPLID, &
                     mpigrp_CPLID, type_grid, glc(1)%cplcompid, idintx)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error in migrating glc mesh for map glc 2 lnd'
                   call shr_sys_abort(subname//' ERROR in migrating glc mesh for map glc 2 lnd')
                endif
                ierr = iMOAB_ComputeCommGraph( mbgxid, mbintxgl, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, &
                     type_grid, type_grid, glc(1)%cplcompid, idintx)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error in computing comm graph for second hop, GLC-LND'
                   call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, GLC-LND')
                endif
                ! prep_lnd_init did not run its glc branch, so the staging and scratch
                ! tags map_glc2lnd_ec writes into do not exist yet; define them
                ! here with the same shared routine prep_lnd_init uses
                call map_glc2lnd_ec_init(mbgxid, mblxid)
             else
                ! reuse the map app and weights loaded by prep_lnd_init
                mapper_Fg2l%src_mbid = mbgxid
                mapper_Fg2l%tgt_mbid = mblxid
                mapper_Fg2l%intx_mbid = mbintxgl
                mapper_Fg2l%src_context = glc(1)%cplcompid
                mapper_Fg2l%intx_context = 100*glc(1)%cplcompid + lnd(1)%cplcompid
                mapper_Fg2l%weight_identifier = wgtIdFg2l
                mapper_Fg2l%mbname = 'mapper_Fg2l'
             endif

          endif ! mblxid and mbgxid

       end if
       call shr_sys_flush(logunit)

    end if

    ! setup needed for either kind of ocn2glc coupling
    if (glc_present .and. (ocn_c2_glctf .or. ocn_c2_glcshelf)) then
       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       o2x_ox => component_get_c2x_cx(ocn(1))
       lsize_o = mct_aVect_lsize(o2x_ox)

       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)

       allocate(o2x_gx(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_gx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_g)
          call mct_aVect_zero(o2x_gx(eoi))
       enddo

       allocate(x2gacc_gx(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(x2gacc_gx(egi), x2g_gx, lsize_g)
          call mct_aVect_zero(x2gacc_gx(egi))
       end do

       x2gacc_gx_cnt = 0
       samegrid_go = .true.
       if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_go = .false.
    end if

    ! setup needed for ocn2glc TF coupling
    ! MOABTODO: give these mappers moab context when the ocn<->glc coupling is ported
    if (glc_present .and. ocn_c2_glctf) then
       call seq_map_mapinit(mapper_So2g_tf, mpicom_CPLID)
    end if

    ! setup needed for ocn2glcshelf coupling
    if (glc_present .and. ocn_c2_glcshelf) then
       call seq_map_mapinit(mapper_So2g_shelf, mpicom_CPLID)
       call seq_map_mapinit(mapper_Fo2g_shelf, mpicom_CPLID)
       !Initialize module-level arrays associated with compute_melt_fluxes
       allocate(oceanTemperature(lsize_g))
       allocate(oceanSalinity(lsize_g))
       allocate(oceanHeatTransferVelocity(lsize_g))
       allocate(oceanSaltTransferVelocity(lsize_g))
       allocate(interfacePressure(lsize_g))
       allocate(iceTemperature(lsize_g))
       allocate(iceTemperatureDistance(lsize_g))
       allocate(iceFloatingMask(lsize_g))
       allocate(outInterfaceSalinity(lsize_g))
       allocate(outInterfaceTemperature(lsize_g))
       allocate(outFreshwaterFlux(lsize_g))
       allocate(outOceanHeatFlux(lsize_g))
       allocate(outIceHeatFlux(lsize_g))
       ! TODO: Can we allocate these only while used or are we worried about performance hit?
       ! TODO: add deallocates!

    end if

    call shr_sys_flush(logunit)

  end subroutine prep_glc_init

  !================================================================================================

  function prep_glc_do_renormalize_smb(infodata) result(do_renormalize_smb)
    ! Returns a logical saying whether we should do the smb renormalization
    logical :: do_renormalize_smb  ! function return value
    !
    ! Arguments
    type (seq_infodata_type) , intent(in) :: infodata

    ! Local variables
    character(len=cl) :: glc_renormalize_smb  ! namelist option saying whether to do smb renormalization
    logical           :: glc_coupled_fluxes   ! does glc send fluxes to other components?
    logical           :: lnd_prognostic       ! is lnd a prognostic component?

    character(len=*), parameter :: subname = '(prep_glc_do_renormalize_smb)'
    !---------------------------------------------------------------

    call seq_infodata_getdata(infodata, &
         glc_renormalize_smb = glc_renormalize_smb, &
         glc_coupled_fluxes  = glc_coupled_fluxes, &
         lnd_prognostic      = lnd_prognostic)

    select case (glc_renormalize_smb)
    case ('on')
       do_renormalize_smb = .true.
    case ('off')
       do_renormalize_smb = .false.
    case ('on_if_glc_coupled_fluxes')
       if (.not. lnd_prognostic) then
          ! Do not renormalize if we're running glc with dlnd (T compsets): In this case
          ! there is no feedback from glc to lnd, and conservation is not important
          do_renormalize_smb = .false.
       else if (.not. glc_coupled_fluxes) then
          ! Do not renormalize if glc isn't sending fluxes to other components: In this
          ! case conservation is not important
          do_renormalize_smb = .false.
       else
          ! lnd_prognostic is true and glc_coupled_fluxes is true
          do_renormalize_smb = .true.
       end if
    case default
       write(logunit,*) subname,' ERROR: unknown value for glc_renormalize_smb: ', &
            trim(glc_renormalize_smb)
       call shr_sys_abort(subname//' ERROR: unknown value for glc_renormalize_smb')
    end select
  end function prep_glc_do_renormalize_smb

  !================================================================================================

  subroutine prep_glc_accum_ocn(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate glc inputs from ocn
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_avect), pointer :: x2g_gx

    character(*), parameter :: subname = '(prep_glc_accum_ocn)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       x2g_gx => component_get_x2c_cx(glc(egi))
       if (x2gacc_gx_cnt == 0) then
          call mct_avect_copy(x2g_gx, x2gacc_gx(egi))
       else
          call mct_avect_accum(x2g_gx, x2gacc_gx(egi))
       endif
    end do
    x2gacc_gx_cnt = x2gacc_gx_cnt + 1
    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_accum_ocn

  !================================================================================================


  subroutine prep_glc_accum_lnd(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate the per-elevation-class land forcing for glc, by reading the
    ! coupler land mesh tags into a plain array (prep_rof accumulation pattern)
    !
    use iMOAB, only : iMOAB_GetDoubleTagStorage
    use iso_c_binding, only : C_NULL_CHAR
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: ierr, ent_type, arrsize
    character(CXX) :: tagname
    character(*), parameter :: subname = '(prep_glc_accum_lnd)'
    !---------------------------------------------------------------

    if (.not. allocated(l2gacc_lm)) return ! moab accumulator not set up

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    tagname = trim(sharedFieldsLndGlc)//C_NULL_CHAR
    arrsize = nflds_lg * lsize_lm
    ent_type = 1 ! cells
    ierr = iMOAB_GetDoubleTagStorage ( mblxid, tagname, arrsize, ent_type, l2x_lm2 )
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' error in getting per-EC lnd fields for glc accumulation')
    endif
    if (l2gacc_lm_cnt == 0) then
       l2gacc_lm = l2x_lm2
    else
       l2gacc_lm = l2gacc_lm + l2x_lm2
    endif
    l2gacc_lm_cnt = l2gacc_lm_cnt + 1
    call t_drvstopf (trim(timer))

  end subroutine prep_glc_accum_lnd

  !================================================================================================

  subroutine prep_glc_accum_avg(timer, lnd2glc_averaged_now)

    !---------------------------------------------------------------
    ! Description
    ! Finalize the accumulation of the land forcing for glc. The averaged values
    ! stay in the l2gacc_lm array (matching the mct l2gacc_lx attribute vector);
    ! prep_glc_calc_l2x_gx stages them into the land mesh tags only for the
    ! duration of the horizontal map, so the l2x tags keep their instantaneous
    ! values for the coupler history.
    !
    ! Note: the driver-mct version also averages the ocn accumulation (x2gacc);
    ! that part of the moab port is deferred with the rest of the ocn<->glc
    ! coupling, so nothing currently consumes what prep_glc_accum_ocn builds up.
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    logical, intent(inout) :: lnd2glc_averaged_now ! Set to .true. if lnd2glc averages were taken this timestep (otherwise left unchanged)
    !
    ! Local Variables
    real(r8) :: ravg ! averaging factor
    character(*), parameter :: subname = '(prep_glc_accum_avg)'
    !---------------------------------------------------------------

    if (.not. allocated(l2gacc_lm)) return ! moab accumulator not set up

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (l2gacc_lm_cnt > 0) then
       lnd2glc_averaged_now = .true.
    end if
    if (l2gacc_lm_cnt > 1) then
       ravg = 1.0_r8/real(l2gacc_lm_cnt, r8)
       l2gacc_lm = l2gacc_lm * ravg
    end if
    l2gacc_lm_cnt_avg = l2gacc_lm_cnt
    l2gacc_lm_cnt = 0
    call t_drvstopf (trim(timer))

  end subroutine prep_glc_accum_avg

  !================================================================================================

  subroutine prep_glc_mrg_ocn(infodata, fractions_gx, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge glc inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_gx(:)
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer :: egi, eoi, efi
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter  :: subname = '(prep_glc_mrg_ocn)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eoi = mod((egi-1),num_inst_ocn) + 1
       efi = mod((egi-1),num_inst_frc) + 1

       x2g_gx => component_get_x2c_cx(glc(egi))
       call prep_glc_merge_ocn_forcing(o2x_gx(eoi), fractions_gx(efi), x2g_gx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg_ocn

  !================================================================================================

  subroutine prep_glc_merge_ocn_forcing( o2x_g, fractions_g, x2g_g )

    !-----------------------------------------------------------------------
    ! Description
    ! "Merge" ocean forcing for glc input.
    !
    ! State fields are copied directly, meaning that averages are taken just over the
    ! ocean-covered portion of the glc domain.
    !
    ! Flux fields are downweighted by landfrac, which effectively sends a 0 flux from the
    ! non-ocean-covered portion of the glc domain.
    !
    ! Arguments
    type(mct_aVect), intent(inout)  :: o2x_g  ! input
    type(mct_aVect), intent(in)     :: fractions_g
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !-----------------------------------------------------------------------

    integer       :: num_flux_fields
    integer       :: num_state_fields
    integer       :: nflds
    integer       :: i,n
    integer       :: mrgstr_index
    integer       :: index_o2x
    integer       :: index_x2g
    integer       :: index_ofrac
    integer       :: lsize
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    character(*), parameter   :: subname = '(prep_glc_merge_ocn_forcing) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2g_g)

    !num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes_from_ocn))
    num_flux_fields = 0
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_tf_states_from_ocn))

    if (first_time) then
       nflds = num_flux_fields + num_state_fields
       allocate(mrgstr(nflds))
    end if

    mrgstr_index = 1

    do i = 1, num_state_fields
       call seq_flds_getField(field, i, seq_flds_x2g_tf_states_from_ocn)
       index_o2x = mct_aVect_indexRA(o2x_g, trim(field))
       index_x2g = mct_aVect_indexRA(x2g_g, trim(field))

       if (first_time) then
          mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
               ' = o2x%'//trim(field)
       end if

       do n = 1, lsize
          x2g_g%rAttr(index_x2g,n) = o2x_g%rAttr(index_o2x,n)
       end do

       mrgstr_index = mrgstr_index + 1
    enddo

    !index_lfrac = mct_aVect_indexRA(fractions_g,"lfrac")
    !do i = 1, num_flux_fields

    !   call seq_flds_getField(field, i, seq_flds_x2g_fluxes_from_lnd)
    !   index_l2x = mct_aVect_indexRA(l2x_g, trim(field))
    !   index_x2g = mct_aVect_indexRA(x2g_g, trim(field))

    !   if (trim(field) == qice_fieldname) then

    !      if (first_time) then
    !         mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
    !              ' = l2x%'//trim(field)
    !      end if

    !      ! treat qice as if it were a state variable, with a simple copy.
    !      do n = 1, lsize
    !         x2g_g%rAttr(index_x2g,n) = l2x_g%rAttr(index_l2x,n)
    !      end do

    !   else
    !      write(logunit,*) subname,' ERROR: Flux fields other than ', &
    !           qice_fieldname, ' currently are not handled in lnd2glc remapping.'
    !      write(logunit,*) '(Attempt to handle flux field <', trim(field), '>.)'
    !      write(logunit,*) 'Substantial thought is needed to determine how to remap other fluxes'
    !      write(logunit,*) 'in a smooth, conservative manner.'
    !      call shr_sys_abort(subname//&
    !           ' ERROR: Flux fields other than qice currently are not handled in lnd2glc remapping.')
    !   endif  ! qice_fieldname

    !   mrgstr_index = mrgstr_index + 1

    !end do

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,nflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif

    first_time = .false.

  end subroutine prep_glc_merge_ocn_forcing


  !================================================================================================

  subroutine prep_glc_calc_o2x_gx(ocn_c2_glctf, ocn_c2_glcshelf, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_gx

    ! Arguments
    character(len=*), intent(in) :: timer
    logical, intent(in) :: ocn_c2_glctf
    logical, intent(in) :: ocn_c2_glcshelf

    character(*), parameter :: subname = '(prep_glc_calc_o2x_gx)'
    ! Local Variables
    integer eoi
    type(mct_avect), pointer :: o2x_ox

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
      o2x_ox => component_get_c2x_cx(ocn(eoi))
      if (ocn_c2_glctf) then
         call seq_map_map(mapper_So2g_tf, o2x_ox, o2x_gx(eoi), &
                       fldlist=seq_flds_x2g_tf_states_from_ocn,norm=.true.)
      end if
      if (ocn_c2_glcshelf) then
         call seq_map_map(mapper_So2g_shelf, o2x_ox, o2x_gx(eoi), &
                       fldlist=seq_flds_x2g_shelf_states_from_ocn,norm=.true.)
      end if
    enddo

    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_o2x_gx

  !================================================================================================


  !================================================================================================

  subroutine prep_glc_calculate_subshelf_boundary_fluxes

    !---------------------------------------------------------------
    ! Description
    ! On the ice sheet grid, calculate shelf boundary fluxes

    use shr_const_mod , only: SHR_CONST_KAPPA_LAND_ICE

    ! Local Variables

    integer :: gsize, n, egi
    type(mct_aVect), pointer :: o2x_ox ! Ocn export, ocn grid, cpl pes
    type(mct_aVect), pointer :: x2g_gx ! Glc import, glc grid, cpl pes
    type(mct_aVect), pointer :: g2x_gx ! Glc import, glc grid, cpl pes

    integer :: index_x2g_So_blt
    integer :: index_x2g_So_bls
    integer :: index_x2g_So_htv
    integer :: index_x2g_So_stv
    integer :: index_x2g_So_rhoeff
    integer :: index_g2x_Sg_tbot
    integer :: index_g2x_Sg_dztbot
    integer :: index_g2x_Sg_lithop
    integer :: index_g2x_Sg_icemask_floating

    integer :: index_g2x_Sg_blis
    integer :: index_g2x_Sg_blit
    integer :: index_g2x_Fogx_qiceho
    integer :: index_g2x_Fogx_qicelo
    integer :: index_x2g_Fogx_qiceli
    integer :: index_x2g_Fogx_qicehi

    character(*), parameter :: subname = '(prep_glc_calculate_subshelf_boundary_fluxes)'
    !---------------------------------------------------------------

    if (.not.(glc_present)) return

    do egi = 1,num_inst_glc

       o2x_ox => component_get_c2x_cx(ocn(egi))
       g2x_gx => component_get_c2x_cx(glc(egi))
       x2g_gx => component_get_x2c_cx(glc(egi))

       !Remap relevant ocean variables to ice sheet grid.
       !Done here instead of in glc-frequency mapping so it happens within ocean coupling interval.
       ! Also could map o2x_ox->o2x_gx(1) but using x2g_gx as destination allows us to see
       ! these fields on the GLC grid of the coupler history file, which helps with debugging.
       call seq_map_map(mapper_So2g_shelf, o2x_ox, x2g_gx, &
       fldlist=seq_flds_x2g_shelf_states_from_ocn,norm=.true.)

       ! inputs to melt flux calculation
       index_x2g_So_blt =    mct_avect_indexra(x2g_gx,'So_blt',perrwith='quiet')
       index_x2g_So_bls =    mct_avect_indexra(x2g_gx,'So_bls',perrwith='quiet')
       index_x2g_So_htv =    mct_avect_indexra(x2g_gx,'So_htv',perrwith='quiet')
       index_x2g_So_stv =    mct_avect_indexra(x2g_gx,'So_stv',perrwith='quiet')
       index_x2g_So_rhoeff = mct_avect_indexra(x2g_gx,'So_rhoeff',perrwith='quiet')

       index_g2x_Sg_tbot =   mct_avect_indexra(g2x_gx,'Sg_tbot',perrwith='quiet')
       index_g2x_Sg_dztbot = mct_avect_indexra(g2x_gx,'Sg_dztbot',perrwith='quiet')
       index_g2x_Sg_lithop = mct_avect_indexra(g2x_gx,'Sg_lithop',perrwith='quiet')
       index_g2x_Sg_icemask_floating = mct_avect_indexra(g2x_gx,'Sg_icemask_floating',perrwith='quiet')

       ! outputs to melt flux calculation
       index_g2x_Sg_blis = mct_avect_indexra(g2x_gx,'Sg_blis',perrwith='quiet')
       index_g2x_Sg_blit = mct_avect_indexra(g2x_gx,'Sg_blit',perrwith='quiet')
       index_g2x_Fogx_qiceho = mct_avect_indexra(g2x_gx,'Fogx_qiceho',perrwith='quiet')
       index_g2x_Fogx_qicelo = mct_avect_indexra(g2x_gx,'Fogx_qicelo',perrwith='quiet')
       index_x2g_Fogx_qiceli = mct_avect_indexra(x2g_gx,'Fogx_qiceli',perrwith='quiet')
       index_x2g_Fogx_qicehi = mct_avect_indexra(x2g_gx,'Fogx_qicehi',perrwith='quiet')

       gsize = mct_aVect_lsize(g2x_gx)

       do n=1,gsize
          !Extract glc and ocn-sourced coupler fields used as input to compute_melt_fluxes to local arrays...

          ! Fields from the ocean, now on the GLC grid
          oceanTemperature(n) =               x2g_gx%rAttr(index_x2g_So_blt,n)
          oceanSalinity(n) =                  x2g_gx%rAttr(index_x2g_So_bls,n)
          oceanHeatTransferVelocity(n) =      x2g_gx%rAttr(index_x2g_So_htv,n)
          oceanSaltTransferVelocity(n) =      x2g_gx%rAttr(index_x2g_So_stv,n)

          ! Fields from the ice sheet model (still on the GLC grid)
          iceTemperature(n) =                 g2x_gx%rAttr(index_g2x_Sg_tbot,n)
          iceTemperatureDistance(n) =         g2x_gx%rAttr(index_g2x_Sg_dztbot,n)
          interfacePressure(n) =              g2x_gx%rAttr(index_g2x_Sg_lithop,n)
          iceFloatingMask(n) =                g2x_gx%rAttr(index_g2x_Sg_icemask_floating,n)

          !... initialize local compute_melt_fluxes output arrays...
          outInterfaceSalinity(n)     =       0.0_r8
          outInterfaceTemperature(n)  =       0.0_r8
          outFreshwaterFlux(n)  =             0.0_r8
          outOceanHeatFlux(n) =               0.0_r8
          outIceHeatFlux(n) =                 0.0_r8
       end do

       !...calculate fluxes...
       call compute_melt_fluxes(oceanTemperature=oceanTemperature,&
                                oceanSalinity=oceanSalinity,&
                                oceanHeatTransferVelocity=oceanHeatTransferVelocity,&
                                oceanSaltTransferVelocity=oceanSaltTransferVelocity,&
                                interfacePressure=interfacePressure,&
                                iceTemperature=iceTemperature,&
                                iceTemperatureDistance=iceTemperatureDistance, &
                                iceFloatingMask=iceFloatingMask, &
                                outInterfaceSalinity=outInterfaceSalinity,&
                                outInterfaceTemperature=outInterfaceTemperature,&
                                outFreshwaterFlux=outFreshwaterFlux,&
                                outOceanHeatFlux=outOceanHeatFlux,&
                                outIceHeatFlux=outIceHeatFlux,&
                                gsize=gsize)

       !...and assign fluxes to glc and ocn-directed coupler fields
       do n=1,gsize
          !Assign outputs from compute_melt_fluxes back into coupler attributes
          g2x_gx%rAttr(index_g2x_Sg_blis,n) =     outInterfaceSalinity(n)    !to ocean
          g2x_gx%rAttr(index_g2x_Sg_blit,n) =     outInterfaceTemperature(n) !to ocean
          g2x_gx%rAttr(index_g2x_Fogx_qiceho,n) = outOceanHeatFlux(n)        !to ocean
          g2x_gx%rAttr(index_g2x_Fogx_qicelo,n)=  outFreshwaterFlux(n)       !to ocean
          x2g_gx%rAttr(index_x2g_Fogx_qicehi,n) = outIceHeatFlux(n)          !to ice sheet
          x2g_gx%rAttr(index_x2g_Fogx_qiceli,n) = -1.0_r8 * outFreshwaterFlux(n)       !to ice sheet
       end do

       !Note: remap ocean-side outputs back onto ocean grid done in call to prep_ocn_shelf_calc_g2x_ox

    end do ! loop over GLC instances

  end subroutine prep_glc_calculate_subshelf_boundary_fluxes

  !================================================================================================

  subroutine prep_glc_zero_fields()

    !---------------------------------------------------------------
    ! Description
    ! Set glc input tags (x2g fields on the coupler-side glc mesh) to zero
    !
    ! See the driver-mct version of this routine for why zeroing is needed for
    ! exact restart tests.

    use iMOAB, only : iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorage
    use seq_comm_mct, only : mbgxid
    use ISO_C_BINDING, only : C_NULL_CHAR
    use shr_kind_mod, only : CXX => shr_kind_CXX
    use shr_string_mod, only : shr_string_listGetNum

    ! Local Variables
    integer :: ierr, ent_type, arrsize, nxflds, lsize_gm
    integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
    real(r8), allocatable :: tmparray(:)
    character(CXX) :: tagname
    character(*), parameter :: subname = '(prep_glc_zero_fields)'
    !---------------------------------------------------------------

    if (mbgxid < 0) return ! nothing to do if the glc coupler mesh does not exist

    ierr = iMOAB_GetMeshInfo ( mbgxid, nvert, nvise, nbl, nsurf, nvisBC )
    if (ierr .ne. 0) then
       write(logunit,*) subname,' cant get size of glc mesh on coupler'
       call shr_sys_abort(subname//' ERROR in getting size of glc mesh on coupler')
    endif
    lsize_gm = nvise(1)
    ent_type = 1 ! cells

    nxflds = shr_string_listGetNum(seq_flds_x2g_fields)
    arrsize = nxflds * lsize_gm
    allocate (tmparray(arrsize))
    tmparray = 0._r8
    tagname = trim(seq_flds_x2g_fields)//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage(mbgxid, tagname, arrsize, ent_type, tmparray)
    if (ierr .ne. 0) then
       write(logunit,*) subname,' cant zero out x2g tags on glc coupler mesh'
       call shr_sys_abort(subname//' cant zero out x2g tags on glc coupler mesh')
    endif
    deallocate (tmparray)

  end subroutine prep_glc_zero_fields

  !================================================================================================

  subroutine prep_glc_calc_l2x_gx(fractions_lx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Builds the glc forcing from the accumulated land fields (this also does the
    ! merge preparation, which driver-mct splits across two routines):
    ! - one batched horizontal map of all the per-elevation-class land fields to the
    !   glc mesh, weighted by lfrac with normalization (the same semantics as the
    !   per-field maps driver-mct does through map_lnd2glc)
    ! - elevation-class vertical interpolation on arrays fetched from the glc mesh tags
    ! - conservative correction of the qice flux (area/aream pre-adjustment and, if
    !   enabled, global smb renormalization)
    ! - results are written directly into the x2g tag names on the glc mesh, which
    !   makes the merge a plain no-op, matching driver-mct's merge (plain copies)
    !
    use iMOAB, only : iMOAB_GetDoubleTagStorage, iMOAB_SetDoubleTagStorage
    use iso_c_binding, only : C_NULL_CHAR
    use map_lnd2glc_mod, only : get_glc_elevation_classes, map_lnd2glc_vertical_interp
    use shr_string_mod, only : shr_string_listGetNum
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: num_flux_fields
    integer :: num_state_fields
    integer :: field_num
    integer :: ierr, ent_type, arrsize, n, ec, kf, nEC
    character(len=cl) :: fieldname
    real(r8), allocatable :: data_lg(:,:)      ! all mapped per-EC fields on the glc mesh
    real(r8), allocatable :: glc_ice_covered(:)
    real(r8), allocatable :: glc_topo(:)
    real(r8), allocatable :: topo_g_EC(:,:)
    real(r8), allocatable :: data_g_EC(:,:)
    real(r8), allocatable :: data_g_bare(:)
    real(r8), allocatable :: data_g(:)         ! downscaled field on the glc mesh
    real(r8), allocatable :: area_g(:)
    real(r8), allocatable :: aream_g(:)
    integer , allocatable :: glc_elevclass(:)
    character(CXX) :: tagname
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

    if ((mblxid < 0) .or. (mbgxid < 0)) return

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)

    ent_type = 1 ! cells

    ! Horizontal map of all per-EC fields at once, weighted by lfrac and normalized.
    ! The averaged accumulation (l2gacc_lm) is staged into the land mesh tags only
    ! for the duration of the map; the instantaneous l2x values are saved first and
    ! restored right after, so the coupler history still shows instantaneous l2x
    ! fields exactly like the mct driver. The attribute vector arguments are
    ! metadata only.
    arrsize = nflds_lg * lsize_lm
    tagname = trim(sharedFieldsLndGlc)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mblxid, tagname, arrsize, ent_type, l2x_lm2)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR saving instantaneous per-EC lnd fields')
    endif
    ierr = iMOAB_SetDoubleTagStorage(mblxid, tagname, arrsize, ent_type, l2gacc_lm)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR staging averaged per-EC lnd fields')
    endif

    call seq_map_map(mapper_Sl2g, l2gacc_lx(1), l2x_gx(1), &
         fldlist=trim(sharedFieldsLndGlc), norm=.true., &
         avwts_s=fractions_lx(1), avwtsfld_s='lfrac')

    ! restore the instantaneous l2x values on the land mesh
    ierr = iMOAB_SetDoubleTagStorage(mblxid, tagname, arrsize, ent_type, l2x_lm2)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR restoring instantaneous per-EC lnd fields')
    endif

    ! fetch the mapped per-EC fields and the glc state needed for the vertical interpolation
    allocate(data_lg(lsize_gm, nflds_lg))
    data_lg = 0._r8
    tagname = trim(sharedFieldsLndGlc)//C_NULL_CHAR
    arrsize = nflds_lg * lsize_gm
    ierr = iMOAB_GetDoubleTagStorage(mbgxid, tagname, arrsize, ent_type, data_lg)
    if (ierr .ne. 0) then
       call shr_sys_abort(subname//' ERROR getting mapped per-EC fields on the glc mesh')
    endif

    allocate(glc_ice_covered(lsize_gm), glc_topo(lsize_gm), glc_elevclass(lsize_gm))
    tagname = trim(Sg_frac_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mbgxid, tagname, lsize_gm, ent_type, glc_ice_covered)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting '//Sg_frac_field)
    tagname = trim(Sg_topo_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mbgxid, tagname, lsize_gm, ent_type, glc_topo)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting '//Sg_topo_field)

    call get_glc_elevation_classes(glc_ice_covered, glc_topo, glc_elevclass)

    nEC = glc_get_num_elevation_classes()
    allocate(topo_g_EC(lsize_gm, nEC))
    allocate(data_g_EC(lsize_gm, nEC))
    allocate(data_g_bare(lsize_gm))
    allocate(data_g(lsize_gm))
    do ec = 1, nEC
       kf = mct_aVect_indexRA(l2gacc_lx(1), 'Sl_topo'//glc_elevclass_as_string(ec))
       topo_g_EC(:,ec) = data_lg(:,kf)
    end do

    ! area arrays needed for the qice conservation correction
    allocate(area_g(lsize_gm), aream_g(lsize_gm))
    tagname = 'area'//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mbgxid, tagname, lsize_gm, ent_type, area_g)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting area on the glc mesh')
    tagname = 'aream'//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mbgxid, tagname, lsize_gm, ent_type, aream_g)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting aream on the glc mesh')

    num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes_from_lnd))
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_states_from_lnd))

    do field_num = 1, num_flux_fields
       call seq_flds_getField(fieldname, field_num, seq_flds_x2g_fluxes_from_lnd)

       if (trim(fieldname) == qice_fieldname) then

          do ec = 1, nEC
             kf = mct_aVect_indexRA(l2gacc_lx(1), trim(fieldname)//glc_elevclass_as_string(ec))
             data_g_EC(:,ec) = data_lg(:,kf)
          end do
          kf = mct_aVect_indexRA(l2gacc_lx(1), trim(fieldname)//glc_elevclass_as_string(0))
          data_g_bare(:) = data_lg(:,kf)

          call map_lnd2glc_vertical_interp(glc_topo, topo_g_EC, data_g_EC, data_g_bare, &
               glc_elevclass, data_g)

          ! Preemptive adjustment for area differences between the glc model and the
          ! coupler: the drv2mdl flux correction on the component side multiplies by
          ! aream/area, so multiply here by area/aream (see the discussion in
          ! prep_glc_map_qice_conservative_lnd2glc)
          do n = 1, lsize_gm
             if (aream_g(n) > 0.0_r8) then
                data_g(n) = data_g(n) * area_g(n)/aream_g(n)
             else
                data_g(n) = 0.0_r8
             endif
          enddo

          if (smb_renormalize) then
             call prep_glc_renormalize_smb(fractions_lx(1), aream_g, data_g)
          end if

          tagname = trim(qice_fieldname)//C_NULL_CHAR
          ierr = iMOAB_SetDoubleTagStorage(mbgxid, tagname, lsize_gm, ent_type, data_g)
          if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR setting '//qice_fieldname)

       else
          write(logunit,*) subname,' ERROR: Flux fields other than ', &
               qice_fieldname, ' currently are not handled in lnd2glc remapping.'
          call shr_sys_abort(subname//&
               ' ERROR: Flux fields other than qice currently are not handled in lnd2glc remapping.')
       endif   ! qice_fieldname
    end do

    do field_num = 1, num_state_fields
       call seq_flds_getField(fieldname, field_num, seq_flds_x2g_states_from_lnd)

       do ec = 1, nEC
          kf = mct_aVect_indexRA(l2gacc_lx(1), trim(fieldname)//glc_elevclass_as_string(ec))
          data_g_EC(:,ec) = data_lg(:,kf)
       end do
       kf = mct_aVect_indexRA(l2gacc_lx(1), trim(fieldname)//glc_elevclass_as_string(0))
       data_g_bare(:) = data_lg(:,kf)

       call map_lnd2glc_vertical_interp(glc_topo, topo_g_EC, data_g_EC, data_g_bare, &
            glc_elevclass, data_g)

       tagname = trim(fieldname)//C_NULL_CHAR
       ierr = iMOAB_SetDoubleTagStorage(mbgxid, tagname, lsize_gm, ent_type, data_g)
       if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR setting '//trim(fieldname))
    end do

    deallocate(data_lg, glc_ice_covered, glc_topo, glc_elevclass)
    deallocate(topo_g_EC, data_g_EC, data_g_bare, data_g)
    deallocate(area_g, aream_g)

    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  subroutine prep_glc_renormalize_smb(fractions_lx, aream_g, qice_g)

    ! Renormalizes the surface mass balance so that the global integral on the glc
    ! grid equals the one on the land grid. Same algorithm as the driver-mct version,
    ! with the land- and glc-side fields fetched from the coupler mesh tags and the
    ! per-EC land fractions recomputed through map_glc2lnd_ec.

    use iMOAB, only : iMOAB_GetDoubleTagStorage
    use iso_c_binding, only : C_NULL_CHAR
    use map_glc2lnd_mod, only : map_glc2lnd_ec
    use shr_mpi_mod, only : shr_mpi_sum, shr_mpi_bcast

    ! Arguments
    type(mct_aVect) , intent(in)    :: fractions_lx ! fractions on the land grid (metadata only)
    real(r8)        , intent(in)    :: aream_g(:)   ! cell areas on glc grid, for mapping
    real(r8)        , intent(inout) :: qice_g(:)    ! qice data on glc grid
    !
    ! Local Variables
    logical :: iamroot
    integer :: ierr, ent_type, n, ec, kf, nEC

    real(r8), allocatable :: aream_l(:)      ! cell areas on land grid, for mapping
    real(r8), allocatable :: lfrac(:)        ! land fraction (lfrin) on land grid
    real(r8), allocatable :: Sg_icemask_l(:) ! icemask on land grid
    real(r8), allocatable :: Sg_icemask_g(:) ! icemask on glc grid
    real(r8), allocatable :: qice_l(:,:)     ! SMB (Flgl_qice) per EC on land grid
    real(r8), allocatable :: frac_l(:,:)     ! EC fractions (Sg_ice_covered) on land grid

    type(mct_aVect) :: av_dum  ! zero-size av for the seq_map_map interface

    real(r8) :: local_accum_on_land_grid
    real(r8) :: global_accum_on_land_grid
    real(r8) :: local_accum_on_glc_grid
    real(r8) :: global_accum_on_glc_grid

    real(r8) :: local_ablat_on_land_grid
    real(r8) :: global_ablat_on_land_grid
    real(r8) :: local_ablat_on_glc_grid
    real(r8) :: global_ablat_on_glc_grid

    real(r8) :: accum_renorm_factor   ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor   ! ratio between global ablation on the two grids

    real(r8) :: effective_area  ! grid cell area multiplied by min(lfrac,Sg_icemask_l)

    character(CXX) :: tagname
    character(*), parameter :: subname = '(prep_glc_renormalize_smb)'
    !---------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    ent_type = 1 ! cells
    nEC = glc_get_num_elevation_classes()

    ! land areas for the mapping and the land fraction basis; note that for E3SM we
    ! are using lfrin instead of lfrac (see the driver-mct version)
    allocate(aream_l(lsize_lm), lfrac(lsize_lm))
    tagname = 'aream'//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mblxid, tagname, lsize_lm, ent_type, aream_l)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting aream on the land mesh')
    tagname = 'lfrin'//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mblxid, tagname, lsize_lm, ent_type, lfrac)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting lfrin on the land mesh')

    ! Map Sg_icemask from the glc mesh to the land mesh tag of the same name (see
    ! the driver-mct version for why this mapping is redone here); the result is
    ! identical to what prep_lnd_calc_g2x_lx produces from the same static g2x data.
    call mct_aVect_init(av_dum, rList=Sg_icemask_field, lsize=0)
    call seq_map_map(mapper_Fg2l, av_dum, av_dum, fldlist=Sg_icemask_field, norm=.true.)
    call mct_aVect_clean(av_dum)
    allocate(Sg_icemask_l(lsize_lm))
    tagname = trim(Sg_icemask_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mblxid, tagname, lsize_lm, ent_type, Sg_icemask_l)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting Sg_icemask on the land mesh')

    ! Map Sg_ice_covered (per elevation class) from glc to land, and return the
    ! normalized per-EC fractions (see the driver-mct version for the rationale)
    allocate(frac_l(lsize_lm, 0:nEC))
    call map_glc2lnd_ec(mapper_Fg2l, &
         frac_field = Sg_frac_field, &
         topo_field = Sg_topo_field, &
         icemask_field = Sg_icemask_field, &
         extra_fields = ' ', &   ! no extra fields
         frac_l_out = frac_l)

    ! qice per elevation class on the land grid, from the averaged accumulator
    allocate(qice_l(lsize_lm, 0:nEC))
    do ec = 0, nEC
       kf = mct_aVect_indexRA(l2gacc_lx(1), qice_fieldname//glc_elevclass_as_string(ec))
       qice_l(:,ec) = l2gacc_lm(:,kf)
    end do

    ! Sum qice over local land grid cells

    local_accum_on_land_grid = 0.0_r8
    local_ablat_on_land_grid = 0.0_r8

    do n = 1, lsize_lm

       effective_area = min(lfrac(n),Sg_icemask_l(n)) * aream_l(n)

       do ec = 0, nEC

          if (qice_l(n,ec) >= 0.0_r8) then
             local_accum_on_land_grid = local_accum_on_land_grid &
                  + effective_area * frac_l(n,ec) * qice_l(n,ec)
          else
             local_ablat_on_land_grid = local_ablat_on_land_grid &
                  + effective_area * frac_l(n,ec) * qice_l(n,ec)
          endif

       enddo  ! ec

    enddo   ! n

    call shr_mpi_sum(local_accum_on_land_grid, &
         global_accum_on_land_grid, &
         mpicom_CPLID, 'accum_l')

    call shr_mpi_sum(local_ablat_on_land_grid, &
         global_ablat_on_land_grid, &
         mpicom_CPLID, 'ablat_l')

    call shr_mpi_bcast(global_accum_on_land_grid, mpicom_CPLID)
    call shr_mpi_bcast(global_ablat_on_land_grid, mpicom_CPLID)

    ! Sum qice_g over local glc grid cells (see the driver-mct version for the
    ! discussion of the areas used here)

    allocate(Sg_icemask_g(size(qice_g)))
    tagname = trim(Sg_icemask_field)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage(mbgxid, tagname, size(qice_g), ent_type, Sg_icemask_g)
    if (ierr .ne. 0) call shr_sys_abort(subname//' ERROR getting Sg_icemask on the glc mesh')

    local_accum_on_glc_grid = 0.0_r8
    local_ablat_on_glc_grid = 0.0_r8

    do n = 1, size(qice_g)

       if (qice_g(n) >= 0.0_r8) then
          local_accum_on_glc_grid = local_accum_on_glc_grid &
               + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       else
          local_ablat_on_glc_grid = local_ablat_on_glc_grid &
               + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       endif

    enddo   ! n

    call shr_mpi_sum(local_accum_on_glc_grid, &
         global_accum_on_glc_grid, &
         mpicom_CPLID, 'accum_g')

    call shr_mpi_sum(local_ablat_on_glc_grid, &
         global_ablat_on_glc_grid, &
         mpicom_CPLID, 'ablat_g')

    call shr_mpi_bcast(global_accum_on_glc_grid, mpicom_CPLID)
    call shr_mpi_bcast(global_ablat_on_glc_grid, mpicom_CPLID)

    ! Renormalize

    if (global_accum_on_glc_grid > 0.0_r8) then
       accum_renorm_factor = global_accum_on_land_grid / global_accum_on_glc_grid
    else
       accum_renorm_factor = 0.0_r8
    endif

    if (global_ablat_on_glc_grid < 0.0_r8) then  ! negative by definition
       ablat_renorm_factor = global_ablat_on_land_grid / global_ablat_on_glc_grid
    else
       ablat_renorm_factor = 0.0_r8
    endif

    if (iamroot) then
       write(logunit,*) 'moab global_accum_on_land_grid = ', global_accum_on_land_grid
       write(logunit,*) 'moab global_accum_on_glc_grid = ', global_accum_on_glc_grid
       write(logunit,*) 'moab accum_renorm_factor = ', accum_renorm_factor
       write(logunit,*) 'moab ablat_renorm_factor = ', ablat_renorm_factor
    endif

    do n = 1, size(qice_g)
       if (qice_g(n) >= 0.0_r8) then
          qice_g(n) = qice_g(n) * accum_renorm_factor
       else
          qice_g(n) = qice_g(n) * ablat_renorm_factor
       endif
    enddo

    deallocate(aream_l)
    deallocate(lfrac)
    deallocate(Sg_icemask_l)
    deallocate(Sg_icemask_g)
    deallocate(qice_l)
    deallocate(frac_l)

  end subroutine prep_glc_renormalize_smb

  !================================================================================================

  subroutine prep_glc_mrg_lnd(infodata, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge the land forcing into the glc input. The driver-mct merge is a plain copy of the
    ! mapped fields into x2g; in the moab path prep_glc_calc_l2x_gx already
    ! writes the downscaled results directly into the x2g tag names on the glc
    ! mesh, so there is nothing left to do here besides optional debug output.
    !
#ifdef MOABDEBUG
    use iMOAB, only : iMOAB_WriteMesh
    use seq_comm_mct, only : num_moab_exports
    use iso_c_binding, only : C_NULL_CHAR
#endif
    !
    ! Arguments
    type(seq_infodata_type) , intent(in) :: infodata
    character(len=*)        , intent(in) :: timer_mrg
    !
    ! Local Variables
#ifdef MOABDEBUG
    integer :: ierr
    character*32 :: outfile, wopts, lnum
#endif
    character(*), parameter :: subname = '(prep_glc_mrg_lnd)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg), barrier=mpicom_CPLID)
#ifdef MOABDEBUG
    if (mbgxid .ge. 0 ) then !  we are on coupler pes, for sure
       write(lnum,"(I0.2)") num_moab_exports
       outfile = 'GlcCplAftMrg'//trim(lnum)//'.h5m'//C_NULL_CHAR
       wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR
       ierr = iMOAB_WriteMesh(mbgxid, trim(outfile), trim(wopts))
       if (ierr .ne. 0) then
          call shr_sys_abort(subname//' error in writing glc mesh after merge')
       endif
    endif
#endif
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_glc_mrg_lnd

  !================================================================================================

  function prep_glc_get_l2gacc_lm()
    real(r8), pointer :: prep_glc_get_l2gacc_lm(:,:)
    prep_glc_get_l2gacc_lm => l2gacc_lm
  end function prep_glc_get_l2gacc_lm

  function prep_glc_get_l2gacc_lm_cnt()
    integer, pointer :: prep_glc_get_l2gacc_lm_cnt
    prep_glc_get_l2gacc_lm_cnt => l2gacc_lm_cnt
  end function prep_glc_get_l2gacc_lm_cnt

  function prep_glc_get_l2gacc_lm_cnt_avg()
    integer, pointer :: prep_glc_get_l2gacc_lm_cnt_avg
    prep_glc_get_l2gacc_lm_cnt_avg => l2gacc_lm_cnt_avg
  end function prep_glc_get_l2gacc_lm_cnt_avg

  function prep_glc_get_sharedFieldsLndGlc()
    character(CXX) :: prep_glc_get_sharedFieldsLndGlc
    prep_glc_get_sharedFieldsLndGlc = sharedFieldsLndGlc
  end function prep_glc_get_sharedFieldsLndGlc

  !================================================================================================

  function prep_glc_get_l2x_gx()
    type(mct_aVect), pointer :: prep_glc_get_l2x_gx(:)
    prep_glc_get_l2x_gx => l2x_gx(:)
  end function prep_glc_get_l2x_gx

  function prep_glc_get_l2gacc_lx()
    type(mct_aVect), pointer :: prep_glc_get_l2gacc_lx(:)
    prep_glc_get_l2gacc_lx => l2gacc_lx(:)
  end function prep_glc_get_l2gacc_lx

  function prep_glc_get_l2gacc_lx_one_instance(lnd_inst)
    integer, intent(in) :: lnd_inst
    type(mct_aVect), pointer :: prep_glc_get_l2gacc_lx_one_instance
    prep_glc_get_l2gacc_lx_one_instance => l2gacc_lx(lnd_inst)
  end function prep_glc_get_l2gacc_lx_one_instance

  function prep_glc_get_l2gacc_lx_cnt()
    integer, pointer :: prep_glc_get_l2gacc_lx_cnt
    prep_glc_get_l2gacc_lx_cnt => l2gacc_lx_cnt
  end function prep_glc_get_l2gacc_lx_cnt

  function prep_glc_get_l2gacc_lx_cnt_avg()
    integer, pointer :: prep_glc_get_l2gacc_lx_cnt_avg
    prep_glc_get_l2gacc_lx_cnt_avg => l2gacc_lx_cnt_avg
  end function prep_glc_get_l2gacc_lx_cnt_avg

  function prep_glc_get_o2x_gx()
    type(mct_aVect), pointer :: prep_glc_get_o2x_gx(:)
    prep_glc_get_o2x_gx => o2x_gx(:)
  end function prep_glc_get_o2x_gx

  function prep_glc_get_x2gacc_gx()
    type(mct_aVect), pointer :: prep_glc_get_x2gacc_gx(:)
    prep_glc_get_x2gacc_gx => x2gacc_gx(:)
  end function prep_glc_get_x2gacc_gx

  function prep_glc_get_x2gacc_gx_cnt()
    integer, pointer :: prep_glc_get_x2gacc_gx_cnt
    prep_glc_get_x2gacc_gx_cnt => x2gacc_gx_cnt
  end function prep_glc_get_x2gacc_gx_cnt

  function prep_glc_get_mapper_Sl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Sl2g
    prep_glc_get_mapper_Sl2g => mapper_Sl2g
  end function prep_glc_get_mapper_Sl2g

  function prep_glc_get_mapper_Fl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Fl2g
    prep_glc_get_mapper_Fl2g => mapper_Fl2g
  end function prep_glc_get_mapper_Fl2g

  function prep_glc_get_mapper_So2g_shelf()
    type(seq_map), pointer :: prep_glc_get_mapper_So2g_shelf
    prep_glc_get_mapper_So2g_shelf=> mapper_So2g_shelf
  end function prep_glc_get_mapper_So2g_shelf

  function prep_glc_get_mapper_Fo2g_shelf()
    type(seq_map), pointer :: prep_glc_get_mapper_Fo2g_shelf
    prep_glc_get_mapper_Fo2g_shelf=> mapper_Fo2g_shelf
  end function prep_glc_get_mapper_Fo2g_shelf

!***********************************************************************
!
!  routine compute_melt_fluxes
!
!> \brief   Computes ocean and ice melt fluxes, etc.
!> \author  Xylar Asay-Davis
!> \date    3/27/2015
!>  This routine computes melt fluxes (melt rate, temperature fluxes
!>  into the ice and the ocean, and salt flux) as well as the interface
!>  temperature and salinity.  This routine expects an ice temperature
!>  in the bottom layer of ice and ocean temperature and salinity in
!>  the top ocean layer as well as the pressure at the ice/ocean interface.
!>
!>  The ocean heat and salt transfer velocities are determined based on
!>  observations of turbulent mixing rates in the under-ice boundary layer.
!>  They should be the product of the friction velocity and a (possibly
!>  spatially variable) non-dimenional transfer coefficient.
!>
!>  The iceTemperatureDistance is the distance between the location
!>  where the iceTemperature is supplied and the ice-ocean interface,
!>  used to compute a temperature gradient.  The ice thermal conductivity,
!>  SHR_CONST_KAPPA_LAND_ICE, is zero for the freezing solution from Holland and Jenkins
!>  (1999) in which the ice is purely insulating.
!
!-----------------------------------------------------------------------

  subroutine compute_melt_fluxes( &
       oceanTemperature, &
       oceanSalinity, &
       oceanHeatTransferVelocity, &
       oceanSaltTransferVelocity, &
       interfacePressure, &
       iceTemperature, &
       iceTemperatureDistance, &
       iceFloatingMask, &
       outInterfaceSalinity, &
       outInterfaceTemperature, &
       outFreshwaterFlux, &
       outOceanHeatFlux, &
       outIceHeatFlux, &
       gsize)

    use shr_const_mod,     only: SHR_CONST_CPICE,  &
                                 SHR_CONST_CPSW,   &
                                 SHR_CONST_LATICE, &
                                 SHR_CONST_RHOICE, &
                                 SHR_CONST_RHOSW,  &
                                 SHR_CONST_DTF_DP, &
                                 SHR_CONST_DTF_DS, &
                                 SHR_CONST_DTF_DPDS, &
                                 SHR_CONST_TF0,    &
                                 SHR_CONST_KAPPA_LAND_ICE

    !-----------------------------------------------------------------
    !
    ! input variables
    !
    !-----------------------------------------------------------------

    real (kind=r8), dimension(:), intent(in) :: &
         oceanTemperature, &          !< Input: ocean temperature in top layer
         oceanSalinity, &             !< Input: ocean salinity in top layer
         oceanHeatTransferVelocity, & !< Input: ocean heat transfer velocity
         oceanSaltTransferVelocity, & !< Input: ocean salt transfer velocity
         interfacePressure, &         !< Input: pressure at the ice-ocean interface
         iceTemperature, &            !< Input: ice temperature in bottom layer
         iceTemperatureDistance       !< Input: distance to ice temperature from ice-ocean interface
    integer, dimension(:), intent(in) :: &
         iceFloatingMask              !< Input: mask of cells that contain floating ice

    integer, intent(in) :: gsize !< Input: number of values in each array

    !-----------------------------------------------------------------
    !
    ! output variables
    !
    !-----------------------------------------------------------------

    real (kind=r8), dimension(:), intent(out) :: &
         outInterfaceSalinity, &    !< Output: ocean salinity at the interface
         outInterfaceTemperature, & !< Output: ice/ocean temperature at the interface
         outFreshwaterFlux, &   !< Output: ocean thickness flux (melt rate)
         outOceanHeatFlux, & !< Output: the temperature flux into the ocean
         outIceHeatFlux      !< Output: the temperature flux into the ice

    !-----------------------------------------------------------------
    !
    ! local variables
    !
    !-----------------------------------------------------------------

    real (kind=r8) :: T0, transferVelocityRatio, Tlatent, nu, a, b, c, eta, &
                         iceHeatFluxCoeff, iceDeltaT, dTf_dS
    integer :: n
    character(*), parameter :: subname = '(compute_melt_fluxes)'

    real (kind=r8), parameter :: minInterfaceSalinity = 0.001_r8

    real (kind=r8), parameter :: referencePressure = 0.0_r8 ! Using reference pressure of 0

    real (kind=r8) :: pressureOffset

    Tlatent = SHR_CONST_LATICE/SHR_CONST_CPSW
    do n = 1, gsize
       if (iceFloatingMask(n) == 0) cycle ! Only calculate on floating cells

       if (oceanHeatTransferVelocity(n) == 0.0_r8) then
          write(logunit,*) 'compute_melt_fluxes ERROR: oceanHeatTransferVelocity value of 0 causes divide by 0 at index ', n
          call shr_sys_abort('compute_melt_fluxes ERROR: oceanHeatTransferVelocity value of 0 causes divide by 0')
       end if

       iceHeatFluxCoeff = SHR_CONST_RHOICE*SHR_CONST_CPICE*SHR_CONST_KAPPA_LAND_ICE/iceTemperatureDistance(n)
       nu = iceHeatFluxCoeff/(SHR_CONST_RHOSW*SHR_CONST_CPSW*oceanHeatTransferVelocity(n))
       pressureOffset = max(interfacePressure(n) - referencePressure, 0.0_r8)
       T0 = SHR_CONST_TF0 + SHR_CONST_DTF_DP * pressureOffset
            !Note: These two terms for T0 are not needed because we are evaluating at salinity=0:
            !+ SHR_CONST_DTF_DS * oceanSalinity(n) + SHR_CONST_DTF_DPDS * pressureOffset * oceanSalinity(n)
       iceDeltaT = T0 - iceTemperature(n)
       dTf_dS = SHR_CONST_DTF_DS + SHR_CONST_DTF_DPDS * pressureOffset

       transferVelocityRatio = oceanSaltTransferVelocity(n)/oceanHeatTransferVelocity(n)

       a = -1.0_r8 * dTf_dS * (1.0_r8 + nu)
       b = transferVelocityRatio*Tlatent - nu*iceDeltaT + oceanTemperature(n) - T0
       c = -transferVelocityRatio*Tlatent*max(oceanSalinity(n), 0.0_r8)
       ! a is non-negative; c is strictly non-positive so we never get imaginary roots.
       ! Since a can be zero, we need a solution of the quadratic equation for 1/Si instead of Si.
       ! Following: https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
       ! Since a and -c are are non-negative, the term in the square root is also always >= |b|.
       ! In all reasonable cases, b will be strictly positive, since transferVelocityRatio*Tlatent ~ 2 C,
       ! T0 ~ -1.8 C and oceanTemperature should never be able to get below about -3 C
       ! As long as either b or both a and c are greater than zero, the strictly non-negative root is
       outInterfaceSalinity(n) = max(-(2.0_r8*c)/(b + sqrt(b**2 - 4.0_r8*a*c)), minInterfaceSalinity)

       outInterfaceTemperature(n) = dTf_dS*outInterfaceSalinity(n)+T0

       outFreshwaterFlux(n) = SHR_CONST_RHOSW*oceanSaltTransferVelocity(n) &
            * (oceanSalinity(n)/outInterfaceSalinity(n) - 1.0_r8)

       ! According to Jenkins et al. (2001), the temperature fluxes into the ocean are:
       !   1. the advection of meltwater into the top layer (or removal for freezing)
       !   2. the turbulent transfer of heat across the boundary layer, based on the termal driving
       outOceanHeatFlux(n) = SHR_CONST_CPSW*(outFreshwaterFlux(n)*outInterfaceTemperature(n) &
            - SHR_CONST_RHOSW*oceanHeatTransferVelocity(n)*(oceanTemperature(n)-outInterfaceTemperature(n)))

       ! the temperature fluxes into the ice are:
       !   1. the advection of ice at the interface temperature out of the domain due to melting
       !      (or in due to freezing)
       !   2. the diffusion (if any) of heat into the ice, based on temperature difference between
       !      the reference point in the ice (either the surface or the middle of the bottom layer)
       !      and the interface
       outIceHeatFlux(n) = -SHR_CONST_CPICE*outFreshwaterFlux(n)*outInterfaceTemperature(n)

       outIceHeatFlux(n) = outIceHeatFlux(n) &
            - iceHeatFluxCoeff*(iceTemperature(n) - outInterfaceTemperature(n))

    end do

  !--------------------------------------------------------------------
  end subroutine compute_melt_fluxes

end module prep_glc_mod
