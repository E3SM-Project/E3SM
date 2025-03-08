module prep_glc_mod

#include "shr_assert.h"
  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_glc, num_inst_lnd, num_inst_frc, &
                              num_inst_ocn
  use seq_comm_mct    , only: CPLID, GLCID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: component_get_dom_cx
  use component_type_mod, only: glc, lnd, ocn
  use glc_elevclass_mod, only : glc_get_num_elevation_classes, glc_elevclass_as_string
  use glc_elevclass_mod, only : glc_all_elevclass_strings, GLC_ELEVCLASS_STRLEN

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_glc_init
  public :: prep_glc_mrg_lnd
  public :: prep_glc_mrg_ocn

  public :: prep_glc_accum_lnd
  public :: prep_glc_accum_ocn
  public :: prep_glc_accum_avg

  public :: prep_glc_calc_l2x_gx
  public :: prep_glc_calc_o2x_gx

  public :: prep_glc_zero_fields

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
  private :: prep_glc_set_g2x_lx_fields
  private :: prep_glc_merge_lnd_forcing
  private :: prep_glc_map_one_state_field_lnd2glc
  private :: prep_glc_map_qice_conservative_lnd2glc
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

  ! Fields needed in the g2x_lx attribute vector used as part of mapping qice from lnd to glc
  character(len=:), allocatable :: g2x_lx_fields

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

  !================================================================================================

contains

  !================================================================================================

  subroutine prep_glc_init(infodata, lnd_c2_glc, ocn_c2_glctf, ocn_c2_glcshelf)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and mapping variables
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

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sl2g'
          end if
          call seq_map_init_rcfile(mapper_Sl2g, lnd(1), glc(1), &
               'seq_maps.rc', 'lnd2glc_smapname:', 'lnd2glc_smaptype:', samegrid_lg, &
               'mapper_Sl2g initialization', esmf_map_flag)

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fl2g'
          end if
          call seq_map_init_rcfile(mapper_Fl2g, lnd(1), glc(1), &
               'seq_maps.rc', 'lnd2glc_fmapname:', 'lnd2glc_fmaptype:', samegrid_lg, &
               'mapper_Fl2g initialization', esmf_map_flag)

          ! We need to initialize our own Fg2l mapper because in some cases (particularly
          ! TG compsets - dlnd forcing CISM) the system doesn't otherwise create a Fg2l
          ! mapper.
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fg2l'
          end if
          call seq_map_init_rcfile(mapper_Fg2l, glc(1), lnd(1), &
               'seq_maps.rc', 'glc2lnd_fmapname:', 'glc2lnd_fmaptype:', samegrid_lg, &
               'mapper_Fg2l initialization', esmf_map_flag)

          call prep_glc_set_g2x_lx_fields()
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
    if (glc_present .and. ocn_c2_glctf) then
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_So2g_tf'
       end if
       call seq_map_init_rcfile(mapper_So2g_tf, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_tf_smapname:','ocn2glc_tf_smaptype:',samegrid_go, &
       'mapper_So2g_tf initialization',esmf_map_flag)
    end if

    ! setup needed for ocn2glcshelf coupling
    if (glc_present .and. ocn_c2_glcshelf) then
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_So2g_shelf'
       end if
       call seq_map_init_rcfile(mapper_So2g_shelf, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_shelf_smapname:','ocn2glc_shelf_smaptype:',samegrid_go, &
       'mapper_So2g_shelf initialization',esmf_map_flag)
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_Fo2g_shelf'
       end if
       call seq_map_init_rcfile(mapper_Fo2g_shelf, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_shelf_fmapname:','ocn2glc_shelf_fmaptype:',samegrid_go, &
       'mapper_Fo2g_shelf initialization',esmf_map_flag)
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

  subroutine prep_glc_set_g2x_lx_fields()

    !---------------------------------------------------------------
    ! Description
    ! Sets the module-level g2x_lx_fields variable.
    !
    ! This gives the fields needed in the g2x_lx attribute vector used as part of mapping
    ! qice from lnd to glc.
    !
    ! Local Variables
    character(len=GLC_ELEVCLASS_STRLEN), allocatable :: all_elevclass_strings(:)
    character(len=:), allocatable :: frac_fields
    character(len=:), allocatable :: topo_fields
    integer :: strlen

    ! 1 is probably enough, but use 10 to be safe, in case the length of the delimiter
    ! changes
    integer, parameter :: extra_len_for_list_merge = 10

    character(len=*), parameter :: subname = '(prep_glc_set_g2x_lx_fields)'
    !---------------------------------------------------------------

    allocate(all_elevclass_strings(0:glc_get_num_elevation_classes()))
    all_elevclass_strings = glc_all_elevclass_strings(include_zero = .true.)
    frac_fields = shr_string_listFromSuffixes( &
         suffixes = all_elevclass_strings, &
         strBase  = Sg_frac_field)
    ! Sg_topo is not actually needed on the land grid in
    ! prep_glc_map_qice_conservative_lnd2glc, but it is required by the current interface
    ! for map_glc2lnd_ec.
    topo_fields = shr_string_listFromSuffixes( &
         suffixes = all_elevclass_strings, &
         strBase  = Sg_topo_field)

    strlen = len_trim(frac_fields) + len_trim(topo_fields) + extra_len_for_list_merge
    allocate(character(len=strlen) :: g2x_lx_fields)
    call shr_string_listMerge(frac_fields, topo_fields, g2x_lx_fields)

  end subroutine prep_glc_set_g2x_lx_fields


  !================================================================================================

  subroutine prep_glc_accum_lnd(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate glc inputs from lnd
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    type(mct_avect), pointer :: l2x_lx

    character(*), parameter :: subname = '(prep_glc_accum_lnd)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       l2x_lx => component_get_c2x_cx(lnd(eli))
       if (l2gacc_lx_cnt == 0) then
          call mct_avect_copy(l2x_lx, l2gacc_lx(eli))
       else
          call mct_avect_accum(l2x_lx, l2gacc_lx(eli))
       endif
    end do
    l2gacc_lx_cnt = l2gacc_lx_cnt + 1
    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_accum_lnd

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


  subroutine prep_glc_accum_avg(timer, lnd2glc_averaged_now)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of glc inputs
    ! Note: There could be separate accum_avg routines for forcing coming
    ! from each component (LND and OCN), but they can be combined here
    ! by taking advantage of l2gacc_lx_cnt and x2gacc_gx_cnt variables
    ! that will only be greater than 0 if corresponding coupling is enabled.
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    logical, intent(inout) :: lnd2glc_averaged_now ! Set to .true. if lnd2glc averages were taken this timestep (otherwise left unchanged)
    !
    ! Local Variables
    integer :: eli, egi
    type(mct_avect), pointer :: x2g_gx

    character(*), parameter :: subname = '(prep_glc_accum_avg)'
    !---------------------------------------------------------------

    ! Accumulation for LND
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (l2gacc_lx_cnt > 0) then
       lnd2glc_averaged_now = .true.
    end if
    if (l2gacc_lx_cnt > 1) then
       do eli = 1,num_inst_lnd
          call mct_avect_avg(l2gacc_lx(eli), l2gacc_lx_cnt)
       end do
    end if
    l2gacc_lx_cnt_avg = l2gacc_lx_cnt
    l2gacc_lx_cnt = 0

    ! Accumulation for OCN
    if (x2gacc_gx_cnt > 1) then
       do egi = 1,num_inst_glc
          ! temporary formation of average
          call mct_avect_avg(x2gacc_gx(egi), x2gacc_gx_cnt)

          ! ***NOTE***THE FOLLOWING ACTUALLY MODIFIES x2g_gx
          x2g_gx => component_get_x2c_cx(glc(egi))
          call mct_avect_copy(x2gacc_gx(egi), x2g_gx)
       enddo
    end if
    x2gacc_gx_cnt = 0

    call t_drvstopf  (trim(timer))

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

  subroutine prep_glc_mrg_lnd(infodata, fractions_gx, timer_mrg)

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
    integer :: egi, eli, efi
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter  :: subname = '(prep_glc_mrg_lnd)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       efi = mod((egi-1),num_inst_frc) + 1

       x2g_gx => component_get_x2c_cx(glc(egi))
       call prep_glc_merge_lnd_forcing(l2x_gx(eli), fractions_gx(efi), x2g_gx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg_lnd

  !================================================================================================

  subroutine prep_glc_merge_lnd_forcing( l2x_g, fractions_g, x2g_g )

    !-----------------------------------------------------------------------
    ! Description
    ! "Merge" land forcing for glc input.
    !
    ! State fields are copied directly, meaning that averages are taken just over the
    ! land-covered portion of the glc domain.
    !
    ! Flux fields are downweighted by landfrac, which effectively sends a 0 flux from the
    ! non-land-covered portion of the glc domain.
    !
    ! Arguments
    type(mct_aVect), intent(inout)  :: l2x_g  ! input
    type(mct_aVect), intent(in)     :: fractions_g
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !-----------------------------------------------------------------------

    integer       :: num_flux_fields
    integer       :: num_state_fields
    integer       :: nflds
    integer       :: i,n
    integer       :: mrgstr_index
    integer       :: index_l2x
    integer       :: index_x2g
    integer       :: index_lfrac
    integer       :: lsize
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    character(*), parameter   :: subname = '(prep_glc_merge_lnd_forcing) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2g_g)

    num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes_from_lnd))
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_states_from_lnd))

    if (first_time) then
       nflds = num_flux_fields + num_state_fields
       allocate(mrgstr(nflds))
    end if

    mrgstr_index = 1

    do i = 1, num_state_fields
       call seq_flds_getField(field, i, seq_flds_x2g_states_from_lnd)
       index_l2x = mct_aVect_indexRA(l2x_g, trim(field))
       index_x2g = mct_aVect_indexRA(x2g_g, trim(field))

       if (first_time) then
          mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
               ' = l2x%'//trim(field)
       end if

       do n = 1, lsize
          x2g_g%rAttr(index_x2g,n) = l2x_g%rAttr(index_l2x,n)
       end do

       mrgstr_index = mrgstr_index + 1
    enddo

    index_lfrac = mct_aVect_indexRA(fractions_g,"lfrac")
    do i = 1, num_flux_fields

       call seq_flds_getField(field, i, seq_flds_x2g_fluxes_from_lnd)
       index_l2x = mct_aVect_indexRA(l2x_g, trim(field))
       index_x2g = mct_aVect_indexRA(x2g_g, trim(field))

       if (trim(field) == qice_fieldname) then

          if (first_time) then
             mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
                  ' = l2x%'//trim(field)
          end if

          ! treat qice as if it were a state variable, with a simple copy.
          do n = 1, lsize
             x2g_g%rAttr(index_x2g,n) = l2x_g%rAttr(index_l2x,n)
          end do

       else
          write(logunit,*) subname,' ERROR: Flux fields other than ', &
               qice_fieldname, ' currently are not handled in lnd2glc remapping.'
          write(logunit,*) '(Attempt to handle flux field <', trim(field), '>.)'
          write(logunit,*) 'Substantial thought is needed to determine how to remap other fluxes'
          write(logunit,*) 'in a smooth, conservative manner.'
          call shr_sys_abort(subname//&
               ' ERROR: Flux fields other than qice currently are not handled in lnd2glc remapping.')
       endif  ! qice_fieldname

       mrgstr_index = mrgstr_index + 1

    end do

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

  end subroutine prep_glc_merge_lnd_forcing


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

  subroutine prep_glc_calc_l2x_gx(fractions_lx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_gx (note that l2x_gx is a local module variable)
    ! Also l2x_gx is really the accumulated l2xacc_lx mapped to l2x_gx
    !
    use shr_string_mod, only : shr_string_listGetNum
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi, eli, efi
    integer :: num_flux_fields
    integer :: num_state_fields
    integer :: field_num
    character(len=cl) :: fieldname
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)

    num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes_from_lnd))
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_states_from_lnd))

    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       efi = mod((egi-1),num_inst_frc) + 1

       do field_num = 1, num_flux_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_fluxes_from_lnd)

          if (trim(fieldname) == qice_fieldname) then

             ! Use a bilinear (Sl2g) mapper, as for states.
             ! The Fg2l mapper is needed to map some glc fields to the land grid
             !  for purposes of conservation.
             call prep_glc_map_qice_conservative_lnd2glc(egi=egi, eli=eli, &
                  fractions_lx = fractions_lx(efi), &
                  mapper_Sl2g = mapper_Sl2g, &
                  mapper_Fg2l = mapper_Fg2l)

          else
             write(logunit,*) subname,' ERROR: Flux fields other than ', &
                  qice_fieldname, ' currently are not handled in lnd2glc remapping.'
             write(logunit,*) '(Attempt to handle flux field <', trim(fieldname), '>.)'
             write(logunit,*) 'Substantial thought is needed to determine how to remap other fluxes'
             write(logunit,*) 'in a smooth, conservative manner.'
             call shr_sys_abort(subname//&
                  ' ERROR: Flux fields other than qice currently are not handled in lnd2glc remapping.')
          endif   ! qice_fieldname

       end do

       do field_num = 1, num_state_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_states_from_lnd)
          call prep_glc_map_one_state_field_lnd2glc(egi=egi, eli=eli, &
               fieldname = fieldname, &
               fractions_lx = fractions_lx(efi), &
               mapper = mapper_Sl2g)
       end do

    enddo   ! egi

    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  subroutine prep_glc_map_one_state_field_lnd2glc(egi, eli, fieldname, fractions_lx, mapper)
    ! Maps a single field from the land grid to the glc grid.
    !
    ! This mapping is not conservative, so should only be used for state fields.
    !
    ! NOTE(wjs, 2017-05-10) We used to map each field separately because each field needed
    ! its own vertical gradient calculator. Now that we don't need vertical gradient
    ! calculators, we may be able to change this to map multiple fields at once, at least
    ! for part of map_lnd2glc.

    use map_lnd2glc_mod, only : map_lnd2glc

    ! Arguments
    integer, intent(in) :: egi  ! glc instance index
    integer, intent(in) :: eli  ! lnd instance index
    character(len=*), intent(in) :: fieldname  ! base name of field to map (without elevation class suffix)
    type(mct_aVect) , intent(in) :: fractions_lx  ! fractions on the land grid, for this frac instance
    type(seq_map), intent(inout) :: mapper
    !
    ! Local Variables
    type(mct_avect), pointer :: g2x_gx    ! glc export, glc grid, cpl pes - allocated in driver
    !---------------------------------------------------------------

    g2x_gx => component_get_c2x_cx(glc(egi))

    call map_lnd2glc(l2x_l = l2gacc_lx(eli), &
         landfrac_l = fractions_lx, &
         g2x_g = g2x_gx, &
         fieldname = fieldname, &
         mapper = mapper, &
         l2x_g = l2x_gx(eli))

  end subroutine prep_glc_map_one_state_field_lnd2glc

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
    ! Set glc inputs to zero
    !
    ! This is appropriate during time intervals when we're not sending valid data to glc.
    ! In principle we shouldn't need to zero the fields at these times (instead, glc
    ! should just ignore the fields at these times). However, some tests (like an ERS or
    ! ERI test that stops the final run segment mid-year) can fail if we don't explicitly
    ! zero the fields, because these x2g fields can then differ upon restart.

    ! Local Variables
    integer :: egi
    type(mct_avect), pointer :: x2g_gx
    !---------------------------------------------------------------


    do egi = 1,num_inst_glc
       x2g_gx => component_get_x2c_cx(glc(egi))
       call mct_aVect_zero(x2g_gx)
    end do
  end subroutine prep_glc_zero_fields

  !================================================================================================

  subroutine prep_glc_map_qice_conservative_lnd2glc(egi, eli, fractions_lx, &
       mapper_Sl2g, mapper_Fg2l)

    ! Maps the surface mass balance field (qice) from the land grid to the glc grid.
    !
    ! Use a smooth, non-conservative (bilinear) mapping, followed by a correction for
    ! conservation.
    !
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit

    use map_lnd2glc_mod, only : map_lnd2glc

    ! Arguments
    integer, intent(in) :: egi  ! glc instance index
    integer, intent(in) :: eli  ! lnd instance index
    type(mct_aVect) , intent(in) :: fractions_lx  ! fractions on the land grid, for this frac instance
    type(seq_map), intent(inout) :: mapper_Sl2g   ! state mapper from land to glc grid; non-conservative
    type(seq_map), intent(inout) :: mapper_Fg2l   ! flux mapper from glc to land grid; conservative
    !
    ! Local Variables
    type(mct_aVect), pointer :: g2x_gx   ! glc export, glc grid

    logical :: iamroot

    !Note: The sums in this subroutine use the coupler areas aream_l and aream_g.
    !      The coupler areas can differ from the native areas area_l and area_g.
    !       (For CISM with a polar stereographic projection, area_g can differ from aream_g
    !       by up to ~10%.)
    !      If so, then the calls to subroutine mct_avect_vecmult in component_mod.F90
    !       (just before and after the call to comp_run) should adjust the SMB fluxes
    !       such that in each grid cell, the native value of area*flux is equal to the
    !       coupler value of aream*flux.  This assumes that the SMB field is contained in
    !       seq_fields l2x_fluxes and seq_fields_x2g_fluxes.

    real(r8), dimension(:), allocatable :: aream_g   ! cell areas on glc grid, for mapping
    real(r8), dimension(:), allocatable :: area_g    ! cell areas on glc grid, according to glc model

    type(mct_ggrid), pointer :: dom_g   ! glc grid info

    integer :: lsize_g   ! number of points on glc grid

    integer :: n
    integer :: km, ka

    real(r8), pointer :: qice_g(:)        ! qice data on glc grid

    !---------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (iamroot) then
       write(logunit,*) ' '
       write(logunit,*) 'In prep_glc_map_qice_conservative_lnd2glc'
       write(logunit,*) 'smb_renormalize = ', smb_renormalize
    endif

    ! Get attribute vector needed for mapping and conservation
    g2x_gx => component_get_c2x_cx(glc(egi))

    ! get grid size
    lsize_g = mct_aVect_lsize(l2x_gx(eli))

    ! allocate and fill area arrays on the glc grid
    ! (Note that we get domain information from instance 1, following what's done in
    ! other parts of the coupler.)
    dom_g => component_get_dom_cx(glc(1))

    allocate(aream_g(lsize_g))
    km = mct_aVect_indexRa(dom_g%data, "aream" )
    aream_g(:) = dom_g%data%rAttr(km,:)

    allocate(area_g(lsize_g))
    ka = mct_aVect_indexRa(dom_g%data, "area" )
    area_g(:) = dom_g%data%rAttr(ka,:)

    ! Map the SMB from the land grid to the glc grid, using a non-conservative state mapper.
    call map_lnd2glc(l2x_l = l2gacc_lx(eli), &
         landfrac_l = fractions_lx, &
         g2x_g = g2x_gx, &
         fieldname = qice_fieldname, &
         mapper = mapper_Sl2g, &
         l2x_g = l2x_gx(eli))

    ! Export the remapped SMB to a local array
    allocate(qice_g(lsize_g))
    call mct_aVect_exportRattr(l2x_gx(eli), trim(qice_fieldname), qice_g)

    ! Make a preemptive adjustment to qice_g to account for area differences between CISM and the coupler.
    !    In component_mod.F90, there is a call to mct_avect_vecmult, which multiplies the fluxes
    !     by aream_g/area_g for conservation purposes. Where CISM areas are larger (area_g > aream_g),
    !     the fluxes are reduced, and where CISM areas are smaller, the fluxes are increased.
    !    As a result, an SMB of 1 m/yr in CLM would be converted to an SMB ranging from
    !     ~0.9 to 1.05 m/yr in CISM (with smaller values where CISM areas are larger, and larger
    !     values where CISM areas are smaller).
    !    Here, to keep CISM values close to the CLM values in the corresponding locations,
    !      we anticipate the later correction and multiply qice_g by area_g/aream_g.
    !     Then the later call to mct_avect_vecmult will bring qice back to the original values
    !       obtained from bilinear remapping.
    !    If Flgl_qice were changed to a state (and not included in seq_flds_x2g_fluxes),
    !     then we could skip this adjustment.
    !
    ! Note that we are free to do this or any other adjustments we want to qice at this
    ! point in the remapping, because the conservation correction will ensure that we
    ! still conserve globally despite these adjustments (and smb_renormalize = .false.
    ! should only be used in cases where conservation doesn't matter anyway).

    do n = 1, lsize_g
       if (aream_g(n) > 0.0_r8) then
          qice_g(n) = qice_g(n) * area_g(n)/aream_g(n)
       else
          qice_g(n) = 0.0_r8
       endif
    enddo

    if (smb_renormalize) then
       call prep_glc_renormalize_smb( &
            eli = eli, &
            fractions_lx = fractions_lx, &
            g2x_gx = g2x_gx, &
            mapper_Fg2l = mapper_Fg2l, &
            aream_g = aream_g, &
            qice_g = qice_g)
    end if

    ! Put the adjusted SMB back into l2x_gx.
    !
    ! If we are doing renormalization, then this is the renormalized SMB. Whether or not
    ! we are doing renormalization, this captures the preemptive adjustment to qice_g to
    ! account for area differences between CISM and the coupler.
    call mct_aVect_importRattr(l2x_gx(eli), qice_fieldname, qice_g)

    ! clean up

    deallocate(aream_g)
    deallocate(area_g)
    deallocate(qice_g)

  end subroutine prep_glc_map_qice_conservative_lnd2glc

  !================================================================================================

  subroutine prep_glc_renormalize_smb(eli, fractions_lx, g2x_gx, mapper_Fg2l, aream_g, qice_g)

    ! Renormalizes surface mass balance (smb, here named qice_g) so that the global
    ! integral on the glc grid is equal to the global integral on the land grid.
    !
    ! This is required for conservation - although conservation is only necessary if we
    ! are running with a fully-interactive, two-way-coupled glc.
    !
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit

    use map_glc2lnd_mod, only : map_glc2lnd_ec

    ! Arguments
    integer         , intent(in)    :: eli          ! lnd instance index
    type(mct_aVect) , intent(in)    :: fractions_lx ! fractions on the land grid, for this frac instance
    type(mct_aVect) , intent(in)    :: g2x_gx       ! glc export, glc grid
    type(seq_map)   , intent(inout) :: mapper_Fg2l  ! flux mapper from glc to land grid; conservative
    real(r8)        , intent(in)    :: aream_g(:)   ! cell areas on glc grid, for mapping
    real(r8)        , intent(inout) :: qice_g(:)    ! qice data on glc grid

    !
    ! Local Variables
    integer :: mpicom
    logical :: iamroot

    type(mct_ggrid), pointer :: dom_l                ! land grid info

    integer :: lsize_l                               ! number of points on land grid
    integer :: lsize_g                               ! number of points on glc grid

    real(r8), dimension(:), allocatable :: aream_l   ! cell areas on land grid, for mapping

    real(r8), pointer :: qice_l(:,:)      ! SMB (Flgl_qice) on land grid
    real(r8), pointer :: frac_l(:,:)      ! EC fractions (Sg_ice_covered) on land grid
    real(r8), pointer :: tmp_field_l(:)   ! temporary field on land grid

    ! The following need to be pointers to satisfy the MCT interface
    ! Note: Sg_icemask defines where the ice sheet model can receive a nonzero SMB from the land model.
    real(r8), pointer :: Sg_icemask_g(:)  ! icemask on glc grid
    real(r8), pointer :: Sg_icemask_l(:)  ! icemask on land grid
    real(r8), pointer :: lfrac(:)         ! land fraction on land grid

    type(mct_aVect) :: g2x_lx            ! glc export, lnd grid (not a pointer: created locally)
    type(mct_avect) :: Sg_icemask_l_av   ! temporary attribute vector holding Sg_icemask on the land grid

    integer :: nEC       ! number of elevation classes
    integer :: n
    integer :: ec
    integer :: km

    ! various strings for building field names
    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: qice_field
    character(len=:), allocatable :: frac_field

    ! local and global sums of accumulation and ablation; used to compute renormalization factors

    real(r8) :: local_accum_on_land_grid
    real(r8) :: global_accum_on_land_grid
    real(r8) :: local_accum_on_glc_grid
    real(r8) :: global_accum_on_glc_grid

    real(r8) :: local_ablat_on_land_grid
    real(r8) :: global_ablat_on_land_grid
    real(r8) :: local_ablat_on_glc_grid
    real(r8) :: global_ablat_on_glc_grid

    ! renormalization factors (should be close to 1, e.g. in range 0.95 to 1.05)
    real(r8) :: accum_renorm_factor   ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor   ! ratio between global ablation on the two grids

    real(r8) :: effective_area  ! grid cell area multiplied by min(lfrac,Sg_icemask_l).
    ! This is the area that can contribute SMB to the ice sheet model.


    !---------------------------------------------------------------

    lsize_g = size(qice_g)
    SHR_ASSERT_FL((size(aream_g) == lsize_g), __FILE__, __LINE__)

    call seq_comm_setptrs(CPLID, mpicom=mpicom)
    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize_l = mct_aVect_lsize(l2gacc_lx(eli))

    ! allocate and fill area arrays on the land grid
    ! (Note that we get domain information from instance 1, following what's done in
    ! other parts of the coupler.)
    dom_l => component_get_dom_cx(lnd(1))

    allocate(aream_l(lsize_l))
    km = mct_aVect_indexRa(dom_l%data, "aream" )
    aream_l(:) = dom_l%data%rAttr(km,:)

    ! Export land fractions from fractions_lx to a local array
    ! Note that for E3SM we are using lfrin instead of lfrac
    allocate(lfrac(lsize_l))
    call mct_aVect_exportRattr(fractions_lx, "lfrin", lfrac)

    ! Map Sg_icemask from the glc grid to the land grid.
    ! This may not be necessary, if Sg_icemask_l has already been mapped from Sg_icemask_g.
    ! It is done here for two reasons:
    ! (1) The mapping will *not* have been done if we are running with dlnd (e.g., a TG case).
    ! (2) Because of coupler lags, the current Sg_icemask_l might not be up to date with
    !     Sg_icemask_g. This probably isn't a problem in practice, but doing the mapping
    !     here ensures the mask is up to date.
    !
    ! This mapping uses the same options as the standard glc -> lnd mapping done in
    ! prep_lnd_calc_g2x_lx. If that mapping ever changed (e.g., changing norm to
    ! .false.), then we should change this mapping, too.
    !
    ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but this
    ! requires some more thought
    call mct_aVect_init(Sg_icemask_l_av, rList = Sg_icemask_field, lsize = lsize_l)
    call seq_map_map(mapper = mapper_Fg2l, &
         av_s = g2x_gx, &
         av_d = Sg_icemask_l_av, &
         fldlist = Sg_icemask_field, &
         norm = .true.)

    ! Export Sg_icemask_l from the temporary attribute vector to a local array
    allocate(Sg_icemask_l(lsize_l))
    call mct_aVect_exportRattr(Sg_icemask_l_av, Sg_icemask_field, Sg_icemask_l)

    ! Clean the temporary attribute vector
    call mct_aVect_clean(Sg_icemask_l_av)

    ! Map Sg_ice_covered from the glc grid to the land grid.
    ! This gives the fields Sg_ice_covered00, Sg_ice_covered01, etc. on the land grid.
    ! These fields are needed to integrate the total SMB on the land grid, for conservation purposes.
    ! As above, the mapping may not be necessary, because Sg_ice_covered might already have been mapped.
    ! However, the mapping will not have been done in a TG case with dlnd, and it might not
    ! be up to date because of coupler lags (though the latter probably isn't a problem
    ! in practice).
    !
    ! Note that, for a case with full two-way coupling, we will only conserve if the
    ! actual land cover used over the course of the year matches these currently-remapped
    ! values. This should generally be the case with the current coupling setup.
    !
    ! One could argue that it would be safer (for conservation purposes) if LND sent its
    ! grid cell average SMB values, or if it sent its own notion of the area in each
    ! elevation class for the purpose of creating grid cell average SMB values here. But
    ! these options cause problems if we're not doing full two-way coupling (e.g., in a TG
    ! case with dlnd, or in the common case where GLC is a diagnostic component that
    ! doesn't cause updates in the glacier areas in LND). In these cases without full
    ! two-way coupling, if we use the LND's notion of the area in each elevation class,
    ! then the conservation corrections would end up correcting for discrepancies in
    ! elevation class areas between LND and GLC, rather than just correcting for
    ! discrepancies arising from the remapping of SMB. (And before you get worried: It
    ! doesn't matter that we are not conserving in these cases without full two-way
    ! coupling, because GLC isn't connected with the rest of the system in terms of energy
    ! and mass in these cases. So in these cases, it's okay that the LND integral computed
    ! here differs from the integral that LND itself would compute.)

    ! Create an attribute vector g2x_lx to hold the mapped fields
    call mct_aVect_init(g2x_lx, rList=g2x_lx_fields, lsize=lsize_l)

    ! Map Sg_ice_covered and Sg_topo from glc to land
    call map_glc2lnd_ec( &
         g2x_g = g2x_gx, &
         frac_field = Sg_frac_field, &
         topo_field = Sg_topo_field, &
         icemask_field = Sg_icemask_field, &
         extra_fields = ' ', &   ! no extra fields
         mapper = mapper_Fg2l, &
         g2x_l = g2x_lx)

    ! Export qice and Sg_ice_covered in each elevation class to local arrays.
    ! Note: qice comes from l2gacc_lx; frac comes from g2x_lx.

    nEC = glc_get_num_elevation_classes()

    allocate(qice_l(lsize_l,0:nEC))
    allocate(frac_l(lsize_l,0:nEC))
    allocate(tmp_field_l(lsize_l))

    do ec = 0, nEC
       elevclass_as_string = glc_elevclass_as_string(ec)

       frac_field = Sg_frac_field // elevclass_as_string    ! Sg_ice_covered01, etc.
       call mct_aVect_exportRattr(g2x_lx, trim(frac_field), tmp_field_l)
       frac_l(:,ec) = tmp_field_l(:)

       qice_field = qice_fieldname // elevclass_as_string    ! Flgl_qice01, etc.
       call mct_aVect_exportRattr(l2gacc_lx(eli), trim(qice_field), tmp_field_l)
       qice_l(:,ec) = tmp_field_l(:)

    enddo

    ! clean the temporary attribute vector g2x_lx
    call mct_aVect_clean(g2x_lx)

    ! Sum qice over local land grid cells

    ! initialize qice sum
    local_accum_on_land_grid = 0.0_r8
    local_ablat_on_land_grid = 0.0_r8

    do n = 1, lsize_l

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
         mpicom, 'accum_l')

    call shr_mpi_sum(local_ablat_on_land_grid, &
         global_ablat_on_land_grid, &
         mpicom, 'ablat_l')

    call shr_mpi_bcast(global_accum_on_land_grid, mpicom)
    call shr_mpi_bcast(global_ablat_on_land_grid, mpicom)

    ! Sum qice_g over local glc grid cells.
    ! Note: This sum uses the coupler areas (aream_g), which differ from the native CISM areas.
    !       But since the original qice_g (from bilinear remapping) has been multiplied by
    !        area_g/aream_g above, this calculation is equivalent to multiplying the original qice_g
    !        by the native CISM areas (area_g).
    !       If Flgl_qice were changed to a state (and not included in seq_flds_x2g_fluxes),
    !        then it would be appropriate to use the native CISM areas in this sum.

    ! Export Sg_icemask from g2x_gx to a local array
    allocate(Sg_icemask_g(lsize_g))
    call mct_aVect_exportRattr(g2x_gx, Sg_icemask_field, Sg_icemask_g)

    local_accum_on_glc_grid = 0.0_r8
    local_ablat_on_glc_grid = 0.0_r8

    do n = 1, lsize_g

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
         mpicom, 'accum_g')

    call shr_mpi_sum(local_ablat_on_glc_grid, &
         global_ablat_on_glc_grid, &
         mpicom, 'ablat_g')

    call shr_mpi_bcast(global_accum_on_glc_grid, mpicom)
    call shr_mpi_bcast(global_ablat_on_glc_grid, mpicom)

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
       write(logunit,*) 'global_accum_on_land_grid = ', global_accum_on_land_grid
       write(logunit,*) 'global_accum_on_glc_grid = ', global_accum_on_glc_grid
       write(logunit,*) 'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,*) 'ablat_renorm_factor = ', ablat_renorm_factor
    endif

    do n = 1, lsize_g
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
    deallocate(tmp_field_l)
    deallocate(qice_l)
    deallocate(frac_l)

  end subroutine prep_glc_renormalize_smb

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
