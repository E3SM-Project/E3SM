module prep_rof_mod

#include "shr_assert.h"
  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: cxx => SHR_KIND_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_lnd, num_inst_rof, num_inst_frc
  use seq_comm_mct,     only: CPLID, ROFID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use shr_log_mod     , only: errMsg => shr_log_errMsg
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: rof, lnd
  use prep_lnd_mod, only: prep_lnd_get_mapper_Fr2l
  use map_lnd2rof_irrig_mod, only: map_lnd2rof_irrig

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_rof_init
  public :: prep_rof_mrg

  public :: prep_rof_accum
  public :: prep_rof_accum_avg

  public :: prep_rof_calc_l2r_rx

  public :: prep_rof_get_l2racc_lx
  public :: prep_rof_get_l2racc_lx_cnt
  public :: prep_rof_get_mapper_Fl2r

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_rof_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Fl2r

  ! attribute vectors
  type(mct_aVect), pointer :: l2r_rx(:)

  ! accumulation variables
  type(mct_aVect), pointer :: l2racc_lx(:)   ! lnd export, lnd grid, cpl pes
  integer        , target  :: l2racc_lx_cnt  ! l2racc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID                            ! MPI cpl communicator

  ! field names and lists, for fields that need to be treated specially
  character(len=*), parameter :: irrig_flux_field = 'Flrl_irrig'
  ! fluxes mapped from lnd to rof that don't need any special handling
  character(CXX) :: lnd2rof_normal_fluxes
  ! whether the model is being run with a separate irrigation field
  logical :: have_irrig_field
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_rof_init(infodata, lnd_c2_rof)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: lnd_c2_rof ! .true.  => lnd to rof coupling on
    !
    ! Local Variables
    integer                     :: lsize_r
    integer                     :: lsize_l
    integer                     :: eli, eri
    logical                     :: samegrid_lr   ! samegrid land and rof
    logical                     :: esmf_map_flag ! .true. => use esmf for mapping
    logical                     :: rof_present   ! .true.  => rof is present
    logical                     :: lnd_present   ! .true.  => lnd is present
    logical                     :: iamroot_CPLID ! .true. => CPLID masterproc
    character(CL)               :: lnd_gnam      ! lnd grid
    character(CL)               :: rof_gnam      ! rof grid
    type(mct_aVect) , pointer   :: l2x_lx
    type(mct_aVect) , pointer   :: x2r_rx
    integer                     :: index_irrig
    character(*)    , parameter :: subname = '(prep_rof_init)'
    character(*)    , parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         rof_present=rof_present       , &
         lnd_present=lnd_present       , &
         lnd_gnam=lnd_gnam             , &
         rof_gnam=rof_gnam             )

    allocate(mapper_Fl2r)

    if (rof_present) then
       x2r_rx => component_get_x2c_cx(rof(1))
       index_irrig = mct_aVect_indexRA(x2r_rx, irrig_flux_field, perrWith='quiet')
       if (index_irrig == 0) then
          have_irrig_field = .false.
       else
          have_irrig_field = .true.
       end if
    else
       ! If rof_present is false, have_irrig_field should be irrelevant; we arbitrarily
       ! set it to false in this case.
       have_irrig_field = .false.
    end if

    if (rof_present .and. lnd_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       lsize_r = mct_aVect_lsize(x2r_rx)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       allocate(l2racc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2r_rx, l2racc_lx(eli), lsize=lsize_l)
          call mct_aVect_zero(l2racc_lx(eli))
       end do
       l2racc_lx_cnt = 0

       allocate(l2r_rx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(l2r_rx(eri), rList=seq_flds_x2r_fields, lsize=lsize_r)
          call mct_avect_zero(l2r_rx(eri))
       end do

       samegrid_lr = .true.
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.

       if (lnd_c2_rof) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fl2r'
          end if
          call seq_map_init_rcfile(mapper_Fl2r, lnd(1), rof(1), &
               'seq_maps.rc','lnd2rof_fmapname:','lnd2rof_fmaptype:',samegrid_lr, &
               string='mapper_Fl2r initialization', esmf_map=esmf_map_flag)

          ! We'll map irrigation specially, so exclude this from the list of l2r fields
          ! that are mapped "normally". Note that the following assumes that all
          ! x2r_fluxes are lnd2rof (as opposed to coming from some other component).
          !
          ! (This listDiff works even if have_irrig_field is false.)
          call shr_string_listDiff( &
               list1 = seq_flds_x2r_fluxes, &
               list2 = irrig_flux_field, &
               listout = lnd2rof_normal_fluxes)
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_rof_init

  !================================================================================================

  subroutine prep_rof_accum(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate land input to river component
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    type(mct_aVect), pointer :: l2x_lx
    character(*), parameter  :: subname = '(prep_rof_accum)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       l2x_lx => component_get_c2x_cx(lnd(eli))
       if (l2racc_lx_cnt == 0) then
          call mct_avect_copy(l2x_lx, l2racc_lx(eli))
       else
          call mct_avect_accum(l2x_lx, l2racc_lx(eli))
       endif
    end do
    l2racc_lx_cnt = l2racc_lx_cnt + 1
    call t_drvstopf (trim(timer))

  end subroutine prep_rof_accum

  !================================================================================================

  subroutine prep_rof_accum_avg(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of land input to river component
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri, eli
    character(*), parameter :: subname = '(prep_rof_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eli = mod((eri-1),num_inst_lnd) + 1
       call mct_avect_avg(l2racc_lx(eli),l2racc_lx_cnt)
    end do
    l2racc_lx_cnt = 0
    call t_drvstopf (trim(timer))

  end subroutine prep_rof_accum_avg

  !================================================================================================

  subroutine prep_rof_mrg(infodata, fractions_rx, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge rof inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_rx(:)
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eri, efi
    type(mct_aVect), pointer :: x2r_rx
    character(*), parameter  :: subname = '(prep_rof_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg), barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       efi = mod((eri-1),num_inst_frc) + 1

       x2r_rx => component_get_x2c_cx(rof(eri))  ! This is actually modifying x2r_rx
       call prep_rof_merge(l2r_rx(eri), fractions_rx(efi), x2r_rx)
    end do
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_rof_mrg

  !================================================================================================

  subroutine prep_rof_merge(l2x_r, fractions_r, x2r_r)

    !-----------------------------------------------------------------------
    ! Description
    ! Merge land rof and ice forcing for rof input
    !
    ! Arguments
    type(mct_aVect),intent(in)    :: l2x_r
    type(mct_aVect),intent(in)    :: fractions_r
    type(mct_aVect),intent(inout) :: x2r_r
    !
    ! Local variables
    integer       :: i
    integer, save :: index_l2x_Flrl_rofsur
    integer, save :: index_l2x_Flrl_rofgwl
    integer, save :: index_l2x_Flrl_rofsub
    integer, save :: index_l2x_Flrl_rofdto
    integer, save :: index_l2x_Flrl_rofi
    integer, save :: index_l2x_Flrl_irrig
    integer, save :: index_x2r_Flrl_rofsur
    integer, save :: index_x2r_Flrl_rofgwl
    integer, save :: index_x2r_Flrl_rofsub
    integer, save :: index_x2r_Flrl_rofdto
    integer, save :: index_x2r_Flrl_rofi
    integer, save :: index_x2r_Flrl_irrig
    integer, save :: index_l2x_Flrl_rofl_16O
    integer, save :: index_l2x_Flrl_rofi_16O
    integer, save :: index_x2r_Flrl_rofl_16O
    integer, save :: index_x2r_Flrl_rofi_16O
    integer, save :: index_l2x_Flrl_rofl_18O
    integer, save :: index_l2x_Flrl_rofi_18O
    integer, save :: index_x2r_Flrl_rofl_18O
    integer, save :: index_x2r_Flrl_rofi_18O
    integer, save :: index_l2x_Flrl_rofl_HDO
    integer, save :: index_l2x_Flrl_rofi_HDO
    integer, save :: index_x2r_Flrl_rofl_HDO
    integer, save :: index_x2r_Flrl_rofi_HDO
    integer, save :: index_lfrac
    logical, save :: first_time = .true.
    logical, save :: flds_wiso_rof = .false.
    real(r8)      :: lfrac
    integer       :: nflds,lsize
    logical       :: iamroot
    character(CL) :: field        ! field string
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(*), parameter   :: subname = '(prep_rof_merge) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2r_r)

    if (first_time) then
       nflds = mct_aVect_nRattr(x2r_r)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2r_r)
          mrgstr(i) = subname//'x2r%'//trim(field)//' ='
       enddo

       index_l2x_Flrl_rofsur = mct_aVect_indexRA(l2x_r,'Flrl_rofsur' )
       index_l2x_Flrl_rofgwl = mct_aVect_indexRA(l2x_r,'Flrl_rofgwl' )
       index_l2x_Flrl_rofsub = mct_aVect_indexRA(l2x_r,'Flrl_rofsub' )
       index_l2x_Flrl_rofdto = mct_aVect_indexRA(l2x_r,'Flrl_rofdto' )
       if (have_irrig_field) then
          index_l2x_Flrl_irrig  = mct_aVect_indexRA(l2x_r,'Flrl_irrig' )
       end if
       index_l2x_Flrl_rofi   = mct_aVect_indexRA(l2x_r,'Flrl_rofi' )
       index_x2r_Flrl_rofsur = mct_aVect_indexRA(x2r_r,'Flrl_rofsur' )
       index_x2r_Flrl_rofgwl = mct_aVect_indexRA(x2r_r,'Flrl_rofgwl' )
       index_x2r_Flrl_rofsub = mct_aVect_indexRA(x2r_r,'Flrl_rofsub' )
       index_x2r_Flrl_rofdto = mct_aVect_indexRA(x2r_r,'Flrl_rofdto' )
       index_x2r_Flrl_rofi   = mct_aVect_indexRA(x2r_r,'Flrl_rofi' )
       if (have_irrig_field) then
          index_x2r_Flrl_irrig  = mct_aVect_indexRA(x2r_r,'Flrl_irrig' )
       end if
       index_l2x_Flrl_rofl_16O = mct_aVect_indexRA(l2x_r,'Flrl_rofl_16O', perrWith='quiet' )

       if ( index_l2x_Flrl_rofl_16O /= 0 ) flds_wiso_rof = .true.
       if ( flds_wiso_rof ) then
          index_l2x_Flrl_rofi_16O = mct_aVect_indexRA(l2x_r,'Flrl_rofi_16O' )
          index_x2r_Flrl_rofl_16O = mct_aVect_indexRA(x2r_r,'Flrl_rofl_16O' )
          index_x2r_Flrl_rofi_16O = mct_aVect_indexRA(x2r_r,'Flrl_rofi_16O' )

          index_l2x_Flrl_rofl_18O = mct_aVect_indexRA(l2x_r,'Flrl_rofl_18O' )
          index_l2x_Flrl_rofi_18O = mct_aVect_indexRA(l2x_r,'Flrl_rofi_18O' )
          index_x2r_Flrl_rofl_18O = mct_aVect_indexRA(x2r_r,'Flrl_rofl_18O' )
          index_x2r_Flrl_rofi_18O = mct_aVect_indexRA(x2r_r,'Flrl_rofi_18O' )

          index_l2x_Flrl_rofl_HDO = mct_aVect_indexRA(l2x_r,'Flrl_rofl_HDO' )
          index_l2x_Flrl_rofi_HDO = mct_aVect_indexRA(l2x_r,'Flrl_rofi_HDO' )
          index_x2r_Flrl_rofl_HDO = mct_aVect_indexRA(x2r_r,'Flrl_rofl_HDO' )
          index_x2r_Flrl_rofi_HDO = mct_aVect_indexRA(x2r_r,'Flrl_rofi_HDO' )
       end if
       index_lfrac = mct_aVect_indexRA(fractions_r,"lfrac")

       index_lfrac = mct_aVect_indexRA(fractions_r,"lfrac")

       mrgstr(index_x2r_Flrl_rofsur) = trim(mrgstr(index_x2r_Flrl_rofsur))//' = '// &
          'lfrac*l2x%Flrl_rofsur'
       mrgstr(index_x2r_Flrl_rofgwl) = trim(mrgstr(index_x2r_Flrl_rofgwl))//' = '// &
          'lfrac*l2x%Flrl_rofgwl'
       mrgstr(index_x2r_Flrl_rofsub) = trim(mrgstr(index_x2r_Flrl_rofsub))//' = '// &
          'lfrac*l2x%Flrl_rofsub'
       mrgstr(index_x2r_Flrl_rofdto) = trim(mrgstr(index_x2r_Flrl_rofdto))//' = '// &
          'lfrac*l2x%Flrl_rofdto'
       mrgstr(index_x2r_Flrl_rofi) = trim(mrgstr(index_x2r_Flrl_rofi))//' = '// &
          'lfrac*l2x%Flrl_rofi'
       if (have_irrig_field) then
          mrgstr(index_x2r_Flrl_irrig) = trim(mrgstr(index_x2r_Flrl_irrig))//' = '// &
               'lfrac*l2x%Flrl_irrig'
       end if
       if ( flds_wiso_rof ) then
          mrgstr(index_x2r_Flrl_rofl_16O) = trim(mrgstr(index_x2r_Flrl_rofl_16O))//' = '// &
             'lfrac*l2x%Flrl_rofl_16O'
          mrgstr(index_x2r_Flrl_rofi_16O) = trim(mrgstr(index_x2r_Flrl_rofi_16O))//' = '// &
             'lfrac*l2x%Flrl_rofi_16O'
          mrgstr(index_x2r_Flrl_rofl_18O) = trim(mrgstr(index_x2r_Flrl_rofl_18O))//' = '// &
             'lfrac*l2x%Flrl_rofl_18O'
          mrgstr(index_x2r_Flrl_rofi_18O) = trim(mrgstr(index_x2r_Flrl_rofi_18O))//' = '// &
             'lfrac*l2x%Flrl_rofi_18O'
          mrgstr(index_x2r_Flrl_rofl_HDO) = trim(mrgstr(index_x2r_Flrl_rofl_HDO))//' = '// &
             'lfrac*l2x%Flrl_rofl_HDO'
          mrgstr(index_x2r_Flrl_rofi_HDO) = trim(mrgstr(index_x2r_Flrl_rofi_HDO))//' = '// &
             'lfrac*l2x%Flrl_rofi_HDO'
       end if
    end if

    do i = 1,lsize
       lfrac = fractions_r%rAttr(index_lfrac,i)
       x2r_r%rAttr(index_x2r_Flrl_rofsur,i) = l2x_r%rAttr(index_l2x_Flrl_rofsur,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofgwl,i) = l2x_r%rAttr(index_l2x_Flrl_rofgwl,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofsub,i) = l2x_r%rAttr(index_l2x_Flrl_rofsub,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofdto,i) = l2x_r%rAttr(index_l2x_Flrl_rofdto,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofi,i) = l2x_r%rAttr(index_l2x_Flrl_rofi,i) * lfrac
       if (have_irrig_field) then
          x2r_r%rAttr(index_x2r_Flrl_irrig,i) = l2x_r%rAttr(index_l2x_Flrl_irrig,i) * lfrac
       end if
       if ( flds_wiso_rof ) then
          x2r_r%rAttr(index_x2r_Flrl_rofl_16O,i) = l2x_r%rAttr(index_l2x_Flrl_rofl_16O,i) * lfrac
          x2r_r%rAttr(index_x2r_Flrl_rofi_16O,i) = l2x_r%rAttr(index_l2x_Flrl_rofi_16O,i) * lfrac
          x2r_r%rAttr(index_x2r_Flrl_rofl_18O,i) = l2x_r%rAttr(index_l2x_Flrl_rofl_18O,i) * lfrac
          x2r_r%rAttr(index_x2r_Flrl_rofi_18O,i) = l2x_r%rAttr(index_l2x_Flrl_rofi_18O,i) * lfrac
          x2r_r%rAttr(index_x2r_Flrl_rofl_HDO,i) = l2x_r%rAttr(index_l2x_Flrl_rofl_HDO,i) * lfrac
          x2r_r%rAttr(index_x2r_Flrl_rofi_HDO,i) = l2x_r%rAttr(index_l2x_Flrl_rofi_HDO,i) * lfrac
       end if
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

  end subroutine prep_rof_merge

  !================================================================================================

  subroutine prep_rof_calc_l2r_rx(fractions_lx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2r_rx (note that l2r_rx is a local module variable)
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri, eli, efi
    type(mct_avect), pointer :: r2x_rx
    type(seq_map)  , pointer :: mapper_Fr2l  ! flux mapper for mapping rof -> lnd
    character(*), parameter :: subname = '(prep_rof_calc_l2r_rx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eli = mod((eri-1),num_inst_lnd) + 1
       efi = mod((eri-1),num_inst_frc) + 1

       ! If the options to this seq_map_map call change (e.g., the use of avwts), similar
       ! changes should be made in map_lnd2rof_irrig.
       call seq_map_map(mapper_Fl2r, l2racc_lx(eli), l2r_rx(eri), &
            fldlist=lnd2rof_normal_fluxes, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')

       if (have_irrig_field) then
          r2x_rx => component_get_c2x_cx(rof(eri))
          mapper_Fr2l => prep_lnd_get_mapper_Fr2l()
          call map_lnd2rof_irrig( &
               l2r_l = l2racc_lx(eli), &
               r2x_r = r2x_rx, &
               irrig_flux_field = irrig_flux_field, &
               avwts_s = fractions_lx(efi), &
               avwtsfld_s = 'lfrin', &
               mapper_Fl2r = mapper_Fl2r, &
               mapper_Fr2l = mapper_Fr2l, &
               l2r_r = l2r_rx(eri))
       end if
    end do
    call t_drvstopf  (trim(timer))

  end subroutine prep_rof_calc_l2r_rx

  !================================================================================================

  function prep_rof_get_l2racc_lx()
    type(mct_aVect), pointer :: prep_rof_get_l2racc_lx(:)
    prep_rof_get_l2racc_lx => l2racc_lx(:)
  end function prep_rof_get_l2racc_lx

  function prep_rof_get_l2racc_lx_cnt()
    integer, pointer :: prep_rof_get_l2racc_lx_cnt
    prep_rof_get_l2racc_lx_cnt => l2racc_lx_cnt
  end function prep_rof_get_l2racc_lx_cnt

  function prep_rof_get_mapper_Fl2r()
    type(seq_map), pointer :: prep_rof_get_mapper_Fl2r
    prep_rof_get_mapper_Fl2r => mapper_Fl2r
  end function prep_rof_get_mapper_Fl2r

end module prep_rof_mod
