module prep_glc_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_glc, num_inst_lnd, num_inst_frc
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
  use component_type_mod, only: component_get_dom_cx   !WHL added to get glc domain info
  use component_type_mod, only: glc, lnd

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_glc_init
  public :: prep_glc_mrg

  public :: prep_glc_accum
  public :: prep_glc_accum_avg

  public :: prep_glc_calc_l2x_gx

  public :: prep_glc_zero_fields

  public :: prep_glc_get_l2x_gx
  public :: prep_glc_get_l2gacc_lx
  public :: prep_glc_get_l2gacc_lx_cnt
  public :: prep_glc_get_mapper_Sl2g
  public :: prep_glc_get_mapper_Fl2g

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_merge
  private :: prep_glc_map_one_field_lnd2glc

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sl2g
  type(seq_map), pointer :: mapper_Fl2g

  !WHL - added a mapper
  type(seq_map), pointer :: mapper_Fg2l

  ! attribute vectors
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator

  !WHL - logic to renormalize the SMB for conservation
  !      Applies only when smb_smooth_downscale = .true.
  !      Should be set to true for 2-way coupled runs with evolving ice sheets.
  !      Probably does not need to be true for 1-way coupling.
  logical, parameter :: smb_renormalize = .true.
!!  logical, parameter :: smb_renormalize = .false.

  ! Name of flux field giving surface mass balance
  character(len=*), parameter :: qice_fieldname = 'Flgl_qice'

  !================================================================================================

contains

  !================================================================================================

  subroutine prep_glc_init(infodata, lnd_c2_glc)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and mapping variables
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: lnd_c2_glc ! .true.  => lnd to glc coupling on
    !
    ! Local Variables
    integer                          :: eli, egi
    integer                          :: lsize_l
    integer                          :: lsize_g
    logical                          :: samegrid_lg   ! samegrid land and glc
    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                          :: glc_present   ! .true. => glc is present
    character(CL)                    :: lnd_gnam      ! lnd grid
    character(CL)                    :: glc_gnam      ! glc grid
    type(mct_avect), pointer         :: l2x_lx
    type(mct_avect), pointer         :: x2g_gx
    character(*), parameter          :: subname = '(prep_glc_init)'
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         glc_present=glc_present       , &
         lnd_gnam=lnd_gnam             , &
         glc_gnam=glc_gnam)

    allocate(mapper_Sl2g)
    allocate(mapper_Fl2g)

    !WHL - added a mapper
    allocate(mapper_Fg2l)

    if (glc_present .and. lnd_c2_glc) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)

       allocate(l2x_gx(num_inst_lnd))
       allocate(l2gacc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_init(l2x_gx(eli), rList=seq_flds_x2g_fields, lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))

          call mct_aVect_init(l2gacc_lx(eli), rList=seq_flds_l2x_fields_to_glc, lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       enddo
       l2gacc_lx_cnt = 0

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

          !WHL - added mapper_Fg2l
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fg2l'
          end if
          call seq_map_init_rcfile(mapper_Fg2l, glc(1), lnd(1), &
               'seq_maps.rc', 'glc2lnd_fmapname:', 'glc2lnd_fmaptype:', samegrid_lg, &
               'mapper_Fg2l initialization', esmf_map_flag)

       end if
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_glc_init

  !================================================================================================

  subroutine prep_glc_accum(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate glc inputs
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    type(mct_avect), pointer :: l2x_lx
    character(*), parameter :: subname = '(prep_glc_accum)'
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

  end subroutine prep_glc_accum

  !================================================================================================

  subroutine prep_glc_accum_avg(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of glc inputs
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    character(*), parameter :: subname = '(prep_glc_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (l2gacc_lx_cnt > 1) then
       do eli = 1,num_inst_lnd
          call mct_avect_avg(l2gacc_lx(eli), l2gacc_lx_cnt)
       end do
    end if
    l2gacc_lx_cnt = 0
    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_accum_avg

  !================================================================================================

  subroutine prep_glc_mrg(infodata, fractions_gx, timer_mrg)

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
    character(*), parameter  :: subname = '(prep_glc_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       efi = mod((egi-1),num_inst_frc) + 1

       x2g_gx => component_get_x2c_cx(glc(egi))
       call prep_glc_merge(l2x_gx(eli), fractions_gx(efi), x2g_gx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg

  !================================================================================================

  subroutine prep_glc_merge( l2x_g, fractions_g, x2g_g )

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
    real(r8)      :: lfrac
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    character(*), parameter   :: subname = '(prep_glc_merge) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2g_g)

    num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes))
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_states))

    if (first_time) then
       nflds = mct_aVect_nRattr(x2g_g)
       if (nflds /= (num_flux_fields + num_state_fields)) then
          write(logunit,*) subname,' ERROR: nflds /= num_flux_fields + num_state_fields: ', &
               nflds, num_flux_fields, num_state_fields
          call shr_sys_abort(subname//' ERROR: nflds /= num_flux_fields + num_state_fields')
       end if

       allocate(mrgstr(nflds))
    end if

    mrgstr_index = 1

    do i = 1, num_state_fields
       call seq_flds_getField(field, i, seq_flds_x2g_states)
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

       call seq_flds_getField(field, i, seq_flds_x2g_fluxes)
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

  end subroutine prep_glc_merge

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

    num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes))
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_states))

    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       efi = mod((egi-1),num_inst_frc) + 1

       do field_num = 1, num_flux_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_fluxes)

          if (trim(fieldname) == qice_fieldname) then

             ! Use a bilinear (Sl2g) mapper, as for states.
             ! The Fg2l mapper is needed to map some glc fields to the land grid
             !  for purposes of conservation.
             call prep_glc_map_qice_conservative_lnd2glc(egi=egi, eli=eli, &
                  fieldname = fieldname, &
                  fractions_lx = fractions_lx(efi), &
                  mapper_Sl2g = mapper_Sl2g, &
                  mapper_Fg2l = mapper_Fg2l)

          else
             write(logunit,*) subname,' ERROR: Flux fields other than ', &
                  qice_fieldname, ' currently are not handled in lnd2glc remapping.'
             write(logunit,*) '(Attempt to handle flux field <', trim(field), '>.)'
             write(logunit,*) 'Substantial thought is needed to determine how to remap other fluxes'
             write(logunit,*) 'in a smooth, conservative manner.'
             call shr_sys_abort(subname//&
                  ' ERROR: Flux fields other than qice currently are not handled in lnd2glc remapping.')
          endif   ! qice_fieldname

       end do

       do field_num = 1, num_state_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_states)
          call prep_glc_map_one_field_lnd2glc(egi=egi, eli=eli, &
               fieldname = fieldname, &
               fractions_lx = fractions_lx(efi), &
               mapper = mapper_Sl2g)
       end do

    enddo   ! egi

    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  subroutine prep_glc_map_one_field_lnd2glc(egi, eli, fieldname, fractions_lx, mapper)
    ! Maps a single field from the land grid to the glc grid.
    !
    ! Note that we remap each field separately because each field needs its own
    ! vertical gradient calculator.

    use vertical_gradient_calculator_2nd_order, only : vertical_gradient_calculator_2nd_order_type
    use vertical_gradient_calculator_factory
    use glc_elevclass_mod, only : glc_get_num_elevation_classes, &
         glc_get_elevclass_bounds, glc_all_elevclass_strings
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

    type(vertical_gradient_calculator_2nd_order_type) :: gradient_calculator
    !---------------------------------------------------------------

    g2x_gx => component_get_c2x_cx(glc(egi))

    gradient_calculator = create_vertical_gradient_calculator_2nd_order( &
         attr_vect = l2gacc_lx(eli), &
         fieldname = fieldname, &
         toponame = 'Sl_topo', &
         elevclass_names = glc_all_elevclass_strings(), &
         elevclass_bounds = glc_get_elevclass_bounds())

    call map_lnd2glc(l2x_l = l2gacc_lx(eli), &
         landfrac_l = fractions_lx, &
         g2x_g = g2x_gx, &
         fieldname = fieldname, &
         gradient_calculator = gradient_calculator, &
         mapper = mapper, &
         l2x_g = l2x_gx(eli))

  end subroutine prep_glc_map_one_field_lnd2glc

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

  subroutine prep_glc_map_qice_conservative_lnd2glc(egi, eli, fieldname, fractions_lx, mapper_Sl2g, mapper_Fg2l)

    ! Maps the surface mass balance field (qice) from the land grid to the glc grid.
    ! Use a smooth, non-conservative (bilinear) mapping, followed by a correction for conservation.

    !WHL - Remove these vertical_gradient use statements after testing
    use vertical_gradient_calculator_2nd_order, only : vertical_gradient_calculator_2nd_order_type
    use vertical_gradient_calculator_factory
    use glc_elevclass_mod, only : glc_get_num_elevation_classes, &
         glc_get_elevclass_bounds, glc_all_elevclass_strings, glc_elevclass_as_string

    use map_lnd2glc_mod, only : map_lnd2glc
    use map_glc2lnd_mod, only : map_glc2lnd_ec

    ! Arguments
    integer, intent(in) :: egi  ! glc instance index
    integer, intent(in) :: eli  ! lnd instance index
    character(len=*), intent(in) :: fieldname  ! base name of field to map (without elevation class suffix)
    type(mct_aVect) , intent(in) :: fractions_lx  ! fractions on the land grid, for this frac instance
    type(seq_map), intent(inout) :: mapper_Sl2g   ! state mapper from land to glc grid; non-conservative
    type(seq_map), intent(inout) :: mapper_Fg2l   ! flux mapper from glc to land grid; conservative
    !
    ! Local Variables
    type(mct_aVect), pointer :: g2x_gx   ! glc export, glc grid
    type(mct_aVect), pointer :: x2l_lx   ! lnd import, lnd grid
    type(mct_aVect), pointer :: g2x_lx   ! glc export, lnd grid

    !WHL - Remove?
    type(vertical_gradient_calculator_2nd_order_type) :: gradient_calculator

    integer :: mpicom    ! mpi comm

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

    real(r8), dimension(:), allocatable :: aream_l   ! cell areas on land grid, for mapping
    real(r8), dimension(:), allocatable :: aream_g   ! cell areas on glc grid, for mapping
    real(r8), dimension(:), allocatable :: area_g    ! cell areas on glc grid, according to glc model

    type(mct_ggrid), pointer :: dom_l   ! land grid info
    type(mct_ggrid), pointer :: dom_g   ! glc grid info

    integer :: lsize_l   ! number of points on land grid
    integer :: lsize_g   ! number of points on glc grid

    integer :: nEC       ! number of elevation classes

    integer :: n, ec
    integer :: km, ka

    real(r8), pointer :: qice_l(:,:)      ! SMB (Flgl_qice) on land grid
    real(r8), pointer :: topo_l(:,:)
    real(r8), pointer :: frac_l(:,:)      ! EC fractions (Sg_ice_covered) on land grid
    real(r8), pointer :: tmp_field_l(:)   ! temporary field on land grid

    ! various strings for building field names
    character(len=:), allocatable :: g2x_fields_from_glc
    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: delimiter
    character(len=:), allocatable :: qice_field
    character(len=:), allocatable :: frac_field
    character(len=:), allocatable :: topo_field

    character(len=*), parameter :: Sg_frac_field = 'Sg_ice_covered'
    character(len=*), parameter :: Sg_topo_field = 'Sg_topo'
    character(len=*), parameter :: Sg_icemask_field = 'Sg_icemask'

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

    ! The following need to be pointers to satisfy the MCT interface
    ! Note: Sg_icemask defines where the ice sheet model can receive a nonzero SMB from the land model.
    real(r8), pointer :: Sg_icemask_g(:)  ! icemask on glc grid
    real(r8), pointer :: Sg_icemask_l(:)  ! icemask on land grid
    real(r8), pointer :: lfrac(:)         ! land fraction on land grid
    real(r8), pointer :: qice_g(:)        ! qice data on glc grid

    ! temporary attribute vectors
    type(mct_avect) :: Sg_icemask_g_av   ! temporary attribute vector holding Sg_icemask on the glc grid
    type(mct_avect) :: Sg_icemask_l_av   ! temporary attribute vector holding Sg_icemask on the land grid

    real(r8) :: effective_area  ! grid cell area multiplied by min(lfrac,Sg_icemask_l).
                                ! This is the area that can contribute SMB to the ice sheet model.

    !WHL - The remaining variables can be removed after testing

    ! parameters for idealized SMB - just for testing
    real(r8), parameter :: q0 = 0.30_r8   ! positive SMB at high elevation
!!    real(r8), parameter :: q0 = 1.0_r8   ! positive SMB at high elevation
    real(r8), parameter :: q1 = 5.e-4     ! SMB gradient with elevation (yr^-1)
!!    real(r8), parameter :: q1 = 0.0_r8     ! SMB gradient with elevation (yr^-1)
    real(r8), parameter :: h0 = 2000._r8  ! elevation above which SMB is constant

!!    logical :: ideal_smb = .true.
    logical :: ideal_smb = .false.
    !---------------------------------------------------------------

    call seq_comm_setptrs(CPLID, mpicom=mpicom)
    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (iamroot) then
       write(logunit,*) ' '
       write(logunit,*) 'In prep_glc_map_qice_conservative_lnd2glc, fieldname =', trim(fieldname)
    endif

    ! Get some attribute vectors needed for mapping and conservation

    g2x_gx => component_get_c2x_cx(glc(egi))
    x2l_lx => component_get_x2c_cx(lnd(eli))

    !WHL - Here is the call to create the vertical gradient calculator.
    !      By default, this is now replaced with linear glint-style interpolation.
    !      Remove after testing.
    gradient_calculator = create_vertical_gradient_calculator_2nd_order( &
         attr_vect = l2gacc_lx(eli), &
         fieldname = fieldname, &
         toponame = 'Sl_topo', &
         elevclass_names = glc_all_elevclass_strings(), &
         elevclass_bounds = glc_get_elevclass_bounds())

    ! get grid sizes
    lsize_l = mct_aVect_lsize(l2gacc_lx(eli))
    lsize_g = mct_aVect_lsize(l2x_gx(eli))

    ! allocate and fill area arrays on the land grid
    dom_l => component_get_dom_cx(lnd(eli))   !WHL - Is eli correct? 

    allocate(aream_l(lsize_l))
    km = mct_aVect_indexRa(dom_l%data, "aream" )
    aream_l(:) = dom_l%data%rAttr(km,:)

    ! allocate and fill area arrays on the glc grid
    dom_g => component_get_dom_cx(glc(egi))   !WHL - Is egi correct? 

    allocate(aream_g(lsize_g))
    km = mct_aVect_indexRa(dom_g%data, "aream" )
    aream_g(:) = dom_g%data%rAttr(km,:)

    allocate(area_g(lsize_g))
    ka = mct_aVect_indexRa(dom_g%data, "area" )
    area_g(:) = dom_g%data%rAttr(ka,:)

    ! Export land fractions from fractions_lx to a local array
    allocate(lfrac(lsize_l))
    call mct_aVect_exportRattr(fractions_lx, "lfrac", lfrac)
    
    !WHL -  Map Sg_icemask from the glc grid to the land grid.
    ! This may not be necessary, if Sg_icemask_l has already been mapped from Sg_icemask_g.
    ! It is done here for two reasons:
    ! (1) The mapping will *not* have been done if we are running with dlnd (e.g., a TG case).
    ! (2) Because of coupler lags, the current Sg_icemask_l might not be up to date with Sg_icemask_g.
    !     Doing the mapping here ensures the mask is up to date.

    ! Export Sg_icemask from g2x_gx to a local array
    allocate(Sg_icemask_g(lsize_g))
    call mct_aVect_exportRattr(g2x_gx, "Sg_icemask", Sg_icemask_g)

    ! Make a temporary attribute vector holding Sg_icemask_g
    call mct_aVect_init(Sg_icemask_g_av, rList = Sg_icemask_field, lsize = lsize_g)
    call mct_aVect_importRattr(Sg_icemask_g_av, Sg_icemask_field, Sg_icemask_g)

    ! Make a temporary attribute vector holding Sg_icemask_l
    allocate(Sg_icemask_l(lsize_l))
    Sg_icemask_l(:) = 0.0_r8
    call mct_aVect_init(Sg_icemask_l_av, rList = Sg_icemask_field, lsize = lsize_l)

    ! Map Sg_icemask from the glc grid to the land grid
    ! This mapping uses the same options as the standard glc -> lnd mapping done in
    ! prep_lnd_calc_g2x_lx. If that mapping ever changed (e.g., introducing an avwts_s
    ! argument), then it's *possible* that we'd want this mapping to change, too.
    call seq_map_map(mapper = mapper_Fg2l, &
         av_s = Sg_icemask_g_av, &
         av_d = Sg_icemask_l_av, &
         fldlist = Sg_icemask_field, &
         norm = .true.)   !WHL - Verify that we want norm = .true.

    ! Export Sg_icemask_l from the temporary attribute vector to a local array
    call mct_aVect_exportRattr(Sg_icemask_l_av, Sg_icemask_field, Sg_icemask_l)

    ! Clean the temporary attribute vectors
    call mct_aVect_clean(Sg_icemask_g_av)
    call mct_aVect_clean(Sg_icemask_l_av)

    !WHL - Map Sg_ice_covered from the glc grid to the land grid.
    !      This gives the fields Sg_ice_covered00, Sg_ice_covered01, etc. on the land grid.
    !      These fields are needed to integrate the total SMB on the land grid, for conservation purposes.
    !      As above, the mapping may not be necessary, because Sg_ice_covered might already have been mapped.
    !      However, the mapping will not have been done in a TG case with dlnd, and it might not
    !       be up to date because of coupler lags.
    !WHL - Should the g2x_fields_from_glc be created just once, at initialization?

    ! Make a list of the frac and topo fields for each EC in g2x_lx

    nEC = glc_get_num_elevation_classes()
    g2x_fields_from_glc = ''
    delimiter = ''

    ! frac fields for each EC
    do ec = 0, nEC
       if (ec > 0) delimiter = ':'
       elevclass_as_string = glc_elevclass_as_string(ec)
       frac_field = Sg_frac_field // elevclass_as_string  ! Sg_ice_covered01, etc.
       g2x_fields_from_glc = g2x_fields_from_glc // delimiter // frac_field
    enddo

    ! topo fields for each EC
    do ec = 0, nEC
       elevclass_as_string = glc_elevclass_as_string(ec)
       topo_field = Sg_topo_field // elevclass_as_string  ! Sg_topo01, etc.
       g2x_fields_from_glc = g2x_fields_from_glc // delimiter // topo_field
    enddo

    ! Create an attribute vector g2x_lx to hold the mapped fields

    allocate(g2x_lx)   ! WHL - allocating here, since this is a temporary local AV
    call mct_aVect_init(g2x_lx, rList=g2x_fields_from_glc, lsize=lsize_l)

    ! Map Sg_ice_covered and Sg_topo from glc to land
    ! Sg_topo is not needed in this subroutine (except for diagnostics and testing), 
    !  but is required by the current interface.
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

    !WHL - It would be possible to export qice_l and frac_l in the same EC loop.
    !      But to support an ideal SMB for testing, we need to first get topo, then set the SMB,
    !      and finally import the SMB to qice_l.

    allocate(qice_l(lsize_l,0:nEC))
    allocate(frac_l(lsize_l,0:nEC))
    allocate(topo_l(lsize_l,0:nEC))   !WHL - not needed in general, but used for ideal SMB
    allocate(tmp_field_l(lsize_l))

    do ec = 0, nEC
       elevclass_as_string = glc_elevclass_as_string(ec)

       frac_field = Sg_frac_field // elevclass_as_string    ! Sg_ice_covered01, etc.
       call mct_aVect_exportRattr(g2x_lx, trim(frac_field), tmp_field_l)
       frac_l(:,ec) = tmp_field_l(:)

       !WHL - topo_l is currently used only for ideal SMB
       topo_field = Sg_topo_field // elevclass_as_string    ! Sg_topo01, etc.
       call mct_aVect_exportRattr(g2x_lx, trim(topo_field), tmp_field_l)
       topo_l(:,ec) = tmp_field_l(:)

    enddo

    ! clean the temporary attribute vector g2x_lx
    call mct_aVect_clean(g2x_lx)

    !WHL debug - option to use an ideal SMB
    if (ideal_smb) then

       do ec = 0, nEC
          elevclass_as_string = glc_elevclass_as_string(ec)

          !WHL - debug - Prescribe inception for class 0 if above topo threshold
          if (ec == 0) then

             ! initialize to zero
             qice_l(:,ec) = 0._r8

             ! set to nonzero value above inception topo threshold
             do n = 1, lsize_l
!                if (topo_l(n,ec) > 500.0_r8) then
                   qice_l(n,ec) = 1.0_r8
                   write(logunit,*) 'Bare ice SMB: n, ec, topo, frac', n, ec, topo_l(n,ec), frac_l(n,ec)
!                else
!                   qice_l(n,ec) = 0.0_r8
!                endif
             enddo

          endif

          if (ec > 0) then
             ! assign qice (m/yr) to each topo value
             do n = 1, lsize_l
                if (topo_l(n,ec) > h0) then
                   qice_l(n,ec) = q0
                else
                   qice_l(n,ec) = q0 - q1*(h0 - topo_l(n,ec)) 
                endif
             enddo  ! n
          endif

          ! convert from m/yr to km/m2/s
          qice_l(:,ec) = qice_l(:,ec) * 917._r8 / 31536000._r8 

          ! Import back into aVect
          qice_field = qice_fieldname // elevclass_as_string    ! Flgl_qice01, etc.
          tmp_field_l(:) = qice_l(:,ec)
          call mct_aVect_importRattr(l2gacc_lx(eli), trim(qice_field), tmp_field_l)
          
       enddo  ! ec

    endif  ! ideal_smb

    ! Export qice in each elevation class to local arrays.

    do ec = 0, nEC
       elevclass_as_string = glc_elevclass_as_string(ec)

       qice_field = qice_fieldname // elevclass_as_string    ! Flgl_qice01, etc.
       call mct_aVect_exportRattr(l2gacc_lx(eli), trim(qice_field), tmp_field_l)
       qice_l(:,ec) = tmp_field_l(:)
    enddo

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

    ! Map the SMB from the land grid to the glc grid, using a non-conservative state mapper.
    call map_lnd2glc(l2x_l = l2gacc_lx(eli), &
         landfrac_l = fractions_lx, &
         g2x_g = g2x_gx, &
         fieldname = fieldname, &
         gradient_calculator = gradient_calculator, &   !WHL - gradient calculator can be removed after testing
         mapper = mapper_Sl2g, &
         l2x_g = l2x_gx(eli))

    ! Export the remapped SMB to a local array
    allocate(qice_g(lsize_g))
    call mct_aVect_exportRattr(l2x_gx(eli), trim(fieldname), qice_g)
    
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

    do n = 1, lsize_g
       if (aream_g(n) > 0.0_r8) then
          qice_g(n) = qice_g(n) * area_g(n)/aream_g(n)
       else
          qice_g(n) = 0.0_r8
       endif
    enddo

    ! Sum qice_g over local glc grid cells.
    ! Note: This sum uses the coupler areas (aream_g), which differ from the native CISM areas.
    !       But since the original qice_g (from bilinear remapping) has been multiplied by
    !        area_g/aream_g above, this calculation is equivalent to multiplying the original qice_g
    !        by the native CISM areas (area_g).
    !       If Flgl_qice were changed to a state (and not included in seq_flds_x2g_fluxes),
    !        then it would be appropriate to use the native CISM areas in this sum.

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

    ! Renormalize for conservation

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
       write(logunit,*) 'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,*) 'ablat_renorm_factor = ', ablat_renorm_factor
    endif

    if (smb_renormalize) then

       do n = 1, lsize_g
          if (qice_g(n) >= 0.0_r8) then
             qice_g(n) = qice_g(n) * accum_renorm_factor
          else
             qice_g(n) = qice_g(n) * ablat_renorm_factor
          endif
       enddo

       ! Put the renormalized SMB back into l2x_gx.
       call mct_aVect_importRattr(l2x_gx(eli), qice_fieldname, qice_g)

    endif  ! smb_renormalize

    ! clean up

    deallocate(aream_l)
    deallocate(aream_g)
    deallocate(lfrac)
    deallocate(Sg_icemask_l)
    deallocate(Sg_icemask_g)
    deallocate(tmp_field_l)
    deallocate(qice_l)
    deallocate(frac_l)
    deallocate(qice_g)

    ! The rest are diagnostic only
    deallocate(topo_l)
    deallocate(area_g)

  end subroutine prep_glc_map_qice_conservative_lnd2glc

  !================================================================================================

  function prep_glc_get_l2x_gx()
    type(mct_aVect), pointer :: prep_glc_get_l2x_gx(:)
    prep_glc_get_l2x_gx => l2x_gx(:)
  end function prep_glc_get_l2x_gx

  function prep_glc_get_l2gacc_lx()
    type(mct_aVect), pointer :: prep_glc_get_l2gacc_lx(:)
    prep_glc_get_l2gacc_lx => l2gacc_lx(:)
  end function prep_glc_get_l2gacc_lx

  function prep_glc_get_l2gacc_lx_cnt()
    integer, pointer :: prep_glc_get_l2gacc_lx_cnt
    prep_glc_get_l2gacc_lx_cnt => l2gacc_lx_cnt
  end function prep_glc_get_l2gacc_lx_cnt

  function prep_glc_get_mapper_Sl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Sl2g
    prep_glc_get_mapper_Sl2g => mapper_Sl2g
  end function prep_glc_get_mapper_Sl2g

  function prep_glc_get_mapper_Fl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Fl2g
    prep_glc_get_mapper_Fl2g => mapper_Fl2g
  end function prep_glc_get_mapper_Fl2g

  !WHL - added an Fo2g mapper
!!  function prep_glc_get_mapper_Fo2g()
!!    type(seq_map), pointer :: prep_glc_get_mapper_Fo2g
!!    prep_glc_get_mapper_Fo2g => mapper_Fo2g
!!  end function prep_glc_get_mapper_Fo2g

end module prep_glc_mod
