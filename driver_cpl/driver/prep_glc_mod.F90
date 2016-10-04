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

  ! attribute vectors 
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
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
       
       if (first_time) then
          mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
               ' = lfrac*l2x%'//trim(field)
       end if
       
       do n = 1, lsize
          lfrac = fractions_g%rAttr(index_lfrac,n)
          x2g_g%rAttr(index_x2g,n) = l2x_g%rAttr(index_l2x,n) * lfrac
       end do

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
          call prep_glc_map_one_field_lnd2glc(egi=egi, eli=eli, &
               fieldname = fieldname, &
               fractions_lx = fractions_lx(efi), &
               mapper = mapper_Fl2g)
       end do

       do field_num = 1, num_state_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_states)
          call prep_glc_map_one_field_lnd2glc(egi=egi, eli=eli, &
               fieldname = fieldname, &
               fractions_lx = fractions_lx(efi), &
               mapper = mapper_Sl2g)
       end do
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  subroutine prep_glc_map_one_field_lnd2glc(egi, eli, fieldname, fractions_lx, mapper)
    ! Maps a single field from the land grid to the glc grid.
    !
    ! Note that we remap each field separately because each field needs its own
    ! vertical gradient calculator.

    use vertical_gradient_calculator_2nd_order, only : vertical_gradient_calculator_2nd_order_type
    use glc_elevclass_mod, only : glc_get_num_elevation_classes
    use map_lnd2glc_mod, only : map_lnd2glc

    ! Arguments
    integer, intent(in) :: egi  ! glc instance index
    integer, intent(in) :: eli  ! lnd instance index
    character(len=*), intent(in) :: fieldname  ! base name of field to map (without elevation class suffix)
    type(mct_aVect) , intent(in) :: fractions_lx  ! fractions on the land grid, for this frac instance
    type(seq_map), intent(inout) :: mapper
    !
    ! Local Variables
    type(mct_avect), pointer :: g2x_gx
    type(vertical_gradient_calculator_2nd_order_type) :: gradient_calculator
    !---------------------------------------------------------------

    g2x_gx => component_get_c2x_cx(glc(egi))

    gradient_calculator = vertical_gradient_calculator_2nd_order_type( &
         attr_vect = l2gacc_lx(eli), &
         fieldname = fieldname, &
         toponame = 'Sl_topo', &
         min_elevation_class = 1, &
         max_elevation_class = glc_get_num_elevation_classes())
    call map_lnd2glc(l2x_l = l2gacc_lx(eli), &
         landfrac_l = fractions_lx, &
         g2x_g = g2x_gx, &
         fieldname = fieldname, &
         gradient_calculator = gradient_calculator, &
         mapper = mapper, &
         l2x_g = l2x_gx(eli))

  end subroutine prep_glc_map_one_field_lnd2glc

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

end module prep_glc_mod
