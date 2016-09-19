module prep_rof_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: cxx => SHR_KIND_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_lnd, num_inst_rof, num_inst_frc
  use seq_comm_mct,     only: CPLID, ROFID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: rof, lnd

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
  private :: prep_rof_map_irrig

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

    if (rof_present .and. lnd_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       x2r_rx => component_get_x2c_cx(rof(1))
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
       index_l2x_Flrl_irrig  = mct_aVect_indexRA(l2x_r,'Flrl_irrig' )
       index_l2x_Flrl_rofi   = mct_aVect_indexRA(l2x_r,'Flrl_rofi' )
       index_x2r_Flrl_rofsur = mct_aVect_indexRA(x2r_r,'Flrl_rofsur' )
       index_x2r_Flrl_rofgwl = mct_aVect_indexRA(x2r_r,'Flrl_rofgwl' )
       index_x2r_Flrl_rofsub = mct_aVect_indexRA(x2r_r,'Flrl_rofsub' )
       index_x2r_Flrl_rofdto = mct_aVect_indexRA(x2r_r,'Flrl_rofdto' )
       index_x2r_Flrl_rofi   = mct_aVect_indexRA(x2r_r,'Flrl_rofi' )
       index_x2r_Flrl_irrig  = mct_aVect_indexRA(x2r_r,'Flrl_irrig' )
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
       mrgstr(index_x2r_Flrl_irrig) = trim(mrgstr(index_x2r_Flrl_irrig))//' = '// &
          'lfrac*l2x%Flrl_irrig'
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
       x2r_r%rAttr(index_x2r_Flrl_irrig,i) = l2x_r%rAttr(index_l2x_Flrl_irrig,i) * lfrac
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
    character(*), parameter :: subname = '(prep_rof_calc_l2r_rx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eli = mod((eri-1),num_inst_lnd) + 1
       efi = mod((eri-1),num_inst_frc) + 1
       call seq_map_map(mapper_Fl2r, l2racc_lx(eli), l2r_rx(eri), &
            fldlist=lnd2rof_normal_fluxes, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
       call prep_rof_map_irrig(eri=eri, eli=eli, &
            avwts_s = fractions_lx(efi), &
            avwtsfld_s = 'lfrin')
    end do
    call t_drvstopf  (trim(timer))

  end subroutine prep_rof_calc_l2r_rx

  !================================================================================================

  subroutine prep_rof_map_irrig(eri, eli, avwts_s, avwtsfld_s)
    !---------------------------------------------------------------
    ! Description
    ! Do custom mapping for the irrigation flux
    !
    ! This mapping first converts the flux to a fraction of volr, then maps this
    ! fraction, then converts the fraction back to a flux on the rof grid.
    !
    ! This uses the field given by irrig_flux_field from l2racc_lx(eli), and sets the
    ! field given by irrig_flux_field in l2r_rx(eri).
    !
    ! In addition, it uses the following fields from other attribute vectors:
    ! - r2x_rx (obtained via component_get_c2x_cx(rof(eri))): Flrr_volrmch
    ! - x2l_lx (obtained via component_get_x2c_cx(lnd(eli))): Flrr_volrmch
    !
    ! Arguments
    integer, intent(in) :: eri  ! rof instance index
    integer, intent(in) :: eli  ! lnd instance index
    type(mct_aVect) , intent(in) :: avwts_s    ! attr vect for source weighting
    character(len=*), intent(in) :: avwtsfld_s ! field in avwts_s to use
    !
    ! Local variables
    integer :: i
    integer :: lsize_l  ! number of land points
    integer :: lsize_r  ! number of rof points
    type(mct_avect), pointer :: x2l_lx
    type(mct_avect), pointer :: r2x_rx
    type(mct_avect) :: irrig_l_av  ! temporary attribute vector holding irrigation fluxes on the land grid
    type(mct_avect) :: irrig_r_av  ! temporary attribute vector holding irrigation_fluxes on the rof grid

    ! The following need to be pointers to satisfy the MCT interface:
    real(r8), pointer :: volr_r(:)        ! river volume on the rof grid
    real(r8), pointer :: volr_l(:)        ! river volume on the land grid
    real(r8), pointer :: irrig_flux_l(:)  ! irrigation flux on the land grid [kg m-2 s-1]
    real(r8), pointer :: irrig_flux_r(:)  ! irrigation flux on the rof grid [kg m-2 s-1]
    real(r8), pointer :: irrig_frac_l(:)  ! irrigation as a fraction of volr, land grid
    real(r8), pointer :: irrig_frac_r(:)  ! irrigation as a fraction of volr, rof grid
    real(r8), pointer :: irrig_excess_l(:) ! irrigation where volr <= 0, land grid
    real(r8), pointer :: irrig_excess_r(:) ! irrigation where volr <= 0, rof grid

    character(len=*), parameter :: volr_field         = 'Flrr_volrmch'
    character(len=*), parameter :: irrig_frac_field   = 'Flrl_irrig_frac'
    character(len=*), parameter :: irrig_excess_field = 'Flrl_irrig_excess'
    character(len=*), parameter :: fields_to_remap = &
         irrig_frac_field // ':' // irrig_excess_field
    !---------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Determine attribute vector sizes
    ! ------------------------------------------------------------------------

    lsize_l = mct_aVect_lsize(l2racc_lx(eli))
    lsize_r = mct_aVect_lsize(l2r_rx(eri))

    ! ------------------------------------------------------------------------
    ! Extract the necessary fields from attribute vectors
    ! ------------------------------------------------------------------------

    allocate(irrig_flux_l(lsize_l))
    call mct_aVect_exportRattr(l2racc_lx(eli), irrig_flux_field, irrig_flux_l)

    r2x_rx => component_get_c2x_cx(rof(eri))
    allocate(volr_r(lsize_r))
    call mct_aVect_exportRattr(r2x_rx, volr_field, volr_r)

    ! TODO(wjs, 2016-09-09) Currently we're getting the already-mapped VOLR on the land
    ! grid. I think it would be more robust (to possible changes in driver sequencing) if
    ! we computed volr_l here, by remapping volr_r from rof -> lnd - rather than relying
    ! on the already-mapped volr_l, which may be out-of-sync with the latest volr in ROF,
    ! depending on when rof -> lnd mapping happens relative to this routine.
    x2l_lx => component_get_x2c_cx(lnd(eli))
    allocate(volr_l(lsize_l))
    call mct_aVect_exportRattr(x2l_lx, volr_field, volr_l)

    ! ------------------------------------------------------------------------
    ! Determine irrigation as a fraction of volr
    !
    ! NOTE(wjs, 2016-09-09) Is it really right to call this a fraction of volr? That would
    ! imply that they have the same units. Would it be more accurate to just call this
    ! 'irrigation normalized by volr'? (Then would need to rename some variables; also
    ! search for 'frac' / 'fraction' in comments.)
    ! ------------------------------------------------------------------------

    ! In order to avoid possible divide by 0, as well as to handle non-sensical negative
    ! volr on the land grid, we divide the land's irrigation flux into two separate flux
    ! components: a component where we have positive volr on the land grid (put in
    ! irrig_frac_l, which is mapped using volr normalization) and a component where we
    ! have zero or negative volr on the land grid (put in irrig_excess_l, which is mapped
    ! as a standard flux). We then remap both of these components to the rof grid, and
    ! then finally add the two components to determine the total irrigation flux on the
    ! rof grid.

    allocate(irrig_frac_l(lsize_l))
    allocate(irrig_excess_l(lsize_l))
    do i = 1, lsize_l
       if (volr_l(i) > 0._r8) then
          irrig_frac_l(i)   = irrig_flux_l(i) / volr_l(i)
          irrig_excess_l(i) = 0._r8
       else
          irrig_frac_l(i)   = 0._r8
          irrig_excess_l(i) = irrig_flux_l(i)
       end if
    end do

    ! ------------------------------------------------------------------------
    ! Map irrigation
    ! ------------------------------------------------------------------------

    call mct_aVect_init(irrig_l_av, rList = fields_to_remap, lsize = lsize_l)
    call mct_aVect_importRattr(irrig_l_av, irrig_frac_field, irrig_frac_l)
    call mct_aVect_importRattr(irrig_l_av, irrig_excess_field, irrig_excess_l)
    call mct_aVect_init(irrig_r_av, rList = fields_to_remap, lsize = lsize_r)

    call seq_map_map(mapper = mapper_Fl2r, &
         av_s = irrig_l_av, &
         av_d = irrig_r_av, &
         fldlist = fields_to_remap, &
         norm = .true., &
         avwts_s = avwts_s, &
         avwtsfld_s = avwtsfld_s)

    allocate(irrig_frac_r(lsize_r))
    allocate(irrig_excess_r(lsize_r))
    call mct_aVect_exportRattr(irrig_r_av, irrig_frac_field, irrig_frac_r)
    call mct_aVect_exportRattr(irrig_r_av, irrig_excess_field, irrig_excess_r)

    ! ------------------------------------------------------------------------
    ! Convert to a total irrigation flux on the ROF grid, and put this in the l2r_rx
    ! attribute vector
    ! ------------------------------------------------------------------------

    allocate(irrig_flux_r(lsize_r))
    do i = 1, lsize_r
       irrig_flux_r(i) = (irrig_frac_r(i) * volr_r(i)) + irrig_excess_r(i)
    end do

    call mct_aVect_importRattr(l2r_rx(eri), irrig_flux_field, irrig_flux_r)

    ! ------------------------------------------------------------------------
    ! Clean up
    ! ------------------------------------------------------------------------

    deallocate(volr_r)
    deallocate(volr_l)
    deallocate(irrig_flux_l)
    deallocate(irrig_flux_r)
    deallocate(irrig_frac_l)
    deallocate(irrig_frac_r)
    deallocate(irrig_excess_l)
    deallocate(irrig_excess_r)
    call mct_aVect_clean(irrig_l_av)
    call mct_aVect_clean(irrig_r_av)

  end subroutine prep_rof_map_irrig

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
