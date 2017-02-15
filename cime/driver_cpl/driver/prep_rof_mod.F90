module prep_rof_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8 
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_lnd, num_inst_rof, num_inst_frc
  use seq_comm_mct,     only: num_inst_ocn
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
  use component_type_mod, only: rof, lnd, ocn

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
  public :: prep_rof_calc_o2r_rx
  
  public :: prep_rof_get_l2racc_lx
  public :: prep_rof_get_l2racc_lx_cnt
  public :: prep_rof_get_o2racc_ox
  public :: prep_rof_get_o2racc_ox_cnt
  public :: prep_rof_get_mapper_Fl2r
  public :: prep_rof_get_mapper_Fo2r

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_rof_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Fl2r
  type(seq_map), pointer :: mapper_Fo2r

  ! attribute vectors 
  type(mct_aVect), pointer :: l2r_rx(:)
  type(mct_aVect), pointer :: o2r_rx(:)

  ! accumulation variables
  type(mct_aVect), pointer :: l2racc_lx(:)   ! lnd export, lnd grid, cpl pes
  integer        , target  :: l2racc_lx_cnt  ! l2racc_lx: number of time samples accumulated

  ! accumulation variables
  type(mct_aVect), pointer :: o2racc_ox(:)   ! lnd export, lnd grid, cpl pes
  integer        , target  :: o2racc_ox_cnt  ! l2racc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID                            ! MPI cpl communicator  
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_rof_init(infodata, lnd_c2_rof, ocn_c2_rof)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: lnd_c2_rof ! .true.  => lnd to rof coupling on
    logical                 , intent(in)    :: ocn_c2_rof ! .true.  => ocn to rof coupling on
    !
    ! Local Variables
    integer                     :: lsize_r
    integer                     :: lsize_l
    integer                     :: lsize_o
    integer                     :: eli, eri, eoi
    logical                     :: samegrid_lr   ! samegrid land and rof
    logical                     :: samegrid_or   ! samegrid ocn and rof
    logical                     :: esmf_map_flag ! .true. => use esmf for mapping
    logical                     :: rof_present   ! .true.  => rof is present
    logical                     :: lnd_present   ! .true.  => lnd is present
    logical                     :: ocn_present   ! .true.  => ocn is present
    logical                     :: iamroot_CPLID ! .true. => CPLID masterproc
    character(CL)               :: lnd_gnam      ! lnd grid
    character(CL)               :: rof_gnam      ! rof grid
    character(CL)               :: ocn_gnam      ! ocn grid
    type(mct_aVect) , pointer   :: l2x_lx 
    type(mct_aVect) , pointer   :: o2x_ox 
    type(mct_aVect) , pointer   :: x2r_rx 
    character(*)    , parameter :: subname = '(prep_rof_init)'
    character(*)    , parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         rof_present=rof_present       , &
         lnd_present=lnd_present       , &
         ocn_present=ocn_present       , &
         lnd_gnam=lnd_gnam             , &
         rof_gnam=rof_gnam             , &
         ocn_gnam=ocn_gnam)

    allocate(mapper_Fl2r)
    allocate(mapper_Fo2r)
    
    if (rof_present .and. lnd_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       x2r_rx => component_get_x2c_cx(rof(1)) 
       lsize_r = mct_aVect_lsize(x2r_rx)

       l2x_lx => component_get_c2x_cx(lnd(1)) 
       lsize_l = mct_aVect_lsize(l2x_lx)
       
       o2x_ox => component_get_c2x_cx(ocn(1)) 
       lsize_o = mct_aVect_lsize(o2x_ox) 

       allocate(l2racc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2r_rx, l2racc_lx(eli), lsize=lsize_l)
          call mct_aVect_zero(l2racc_lx(eli))
       end do
       l2racc_lx_cnt = 0
       
       allocate(o2racc_ox(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_initSharedFields(o2x_ox, x2r_rx, o2racc_ox(eoi), lsize=lsize_o)
          call mct_aVect_zero(o2racc_ox(eoi))
       end do
       
       o2racc_ox_cnt = 0

       allocate(l2r_rx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(l2r_rx(eri), rList=seq_flds_x2r_fields, lsize=lsize_r)
          call mct_avect_zero(l2r_rx(eri))
       end do
       
       allocate(o2r_rx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(o2r_rx(eri), rList=seq_flds_x2r_fields, lsize=lsize_r)
          call mct_avect_zero(o2r_rx(eri))
       end do

       samegrid_lr = .true. 
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.
       
       samegrid_or = .true. 
       if (trim(ocn_gnam) /= trim(rof_gnam)) samegrid_or = .false.

       if (lnd_c2_rof) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fl2r'
          end if
          call seq_map_init_rcfile(mapper_Fl2r, lnd(1), rof(1), &
               'seq_maps.rc','lnd2rof_fmapname:','lnd2rof_fmaptype:',samegrid_lr, &
               string='mapper_Fl2r initialization', esmf_map=esmf_map_flag)
       endif

       write(logunit,*) 'samegrid_or=',samegrid_or

       if (ocn_c2_rof) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fo2r'
          end if
          call seq_map_init_rcfile(mapper_Fo2r, ocn(1), rof(1), &
               'seq_maps.rc','ocn2rof_fmapname:','ocn2rof_fmaptype:',samegrid_or, &
               string='mapper_Fo2r initialization', esmf_map=esmf_map_flag)
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
    integer :: eli, eoi
    type(mct_aVect), pointer :: l2x_lx
    type(mct_aVect), pointer :: o2x_ox
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
    
    do eoi = 1,num_inst_ocn
       o2x_ox => component_get_c2x_cx(ocn(eoi)) 
       if (o2racc_ox_cnt == 0) then
          call mct_avect_copy(o2x_ox, o2racc_ox(eoi))
       else
          call mct_avect_accum(o2x_ox, o2racc_ox(eoi))
       endif
    end do
    o2racc_ox_cnt = o2racc_ox_cnt + 1
    
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
    integer :: eri, eli, eoi
    character(*), parameter :: subname = '(prep_rof_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eli = mod((eri-1),num_inst_lnd) + 1
       call mct_avect_avg(l2racc_lx(eli),l2racc_lx_cnt)
    end do
    l2racc_lx_cnt = 0 
    
    do eri = 1,num_inst_rof
       eoi = mod((eri-1),num_inst_ocn) + 1
       call mct_avect_avg(o2racc_ox(eoi),o2racc_ox_cnt)
    end do
    o2racc_ox_cnt = 0
    
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
       call prep_rof_merge(l2r_rx(eri), o2r_rx(eri), fractions_rx(efi), x2r_rx)
    end do
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_rof_mrg

  !================================================================================================

  subroutine prep_rof_merge(l2x_r, o2x_r, fractions_r, x2r_r)

    !----------------------------------------------------------------------- 
    ! Description
    ! Merge land rof and ice forcing for rof input
    !
    ! Arguments
    type(mct_aVect),intent(in)    :: l2x_r
    type(mct_aVect),intent(in)    :: o2x_r
    type(mct_aVect),intent(in)    :: fractions_r
    type(mct_aVect),intent(inout) :: x2r_r
    !
    ! Local variables
    integer       :: i
    integer, save :: index_l2x_Flrl_rofsur
    integer, save :: index_l2x_Flrl_rofgwl
    integer, save :: index_l2x_Flrl_rofsub
    integer, save :: index_l2x_Flrl_rofi
    integer, save :: index_o2x_Floo_qice
    integer, save :: index_x2r_Flrl_rofsur
    integer, save :: index_x2r_Flrl_rofgwl
    integer, save :: index_x2r_Flrl_rofsub
    integer, save :: index_x2r_Flrl_rofi
    integer, save :: index_lfrac
    logical, save :: first_time = .true.
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
       index_l2x_Flrl_rofi   = mct_aVect_indexRA(l2x_r,'Flrl_rofi' )
       index_o2x_Floo_qice   = mct_aVect_indexRA(o2x_r,'Floo_qice' )
       index_x2r_Flrl_rofsur = mct_aVect_indexRA(x2r_r,'Flrl_rofsur' )
       index_x2r_Flrl_rofgwl = mct_aVect_indexRA(x2r_r,'Flrl_rofgwl' )
       index_x2r_Flrl_rofsub = mct_aVect_indexRA(x2r_r,'Flrl_rofsub' )
       index_x2r_Flrl_rofi   = mct_aVect_indexRA(x2r_r,'Flrl_rofi' )
       index_lfrac = mct_aVect_indexRA(fractions_r,"lfrac")

       mrgstr(index_x2r_Flrl_rofsur) = trim(mrgstr(index_x2r_Flrl_rofsur))//' = '// &
          'lfrac*l2x%Flrl_rofsur'
       mrgstr(index_x2r_Flrl_rofgwl) = trim(mrgstr(index_x2r_Flrl_rofgwl))//' = '// &
          'lfrac*l2x%Flrl_rofgwl'
       mrgstr(index_x2r_Flrl_rofsub) = trim(mrgstr(index_x2r_Flrl_rofsub))//' = '// &
          'lfrac*l2x%Flrl_rofsub'
       mrgstr(index_x2r_Flrl_rofi) = trim(mrgstr(index_x2r_Flrl_rofi))//' = '// &
          'lfrac*(l2x%Flrl_rofi - o2x%Floo_qice)'
    end if

    do i = 1,lsize
       lfrac = fractions_r%rAttr(index_lfrac,i)
       x2r_r%rAttr(index_x2r_Flrl_rofsur,i) = l2x_r%rAttr(index_l2x_Flrl_rofsur,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofgwl,i) = l2x_r%rAttr(index_l2x_Flrl_rofgwl,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofsub,i) = l2x_r%rAttr(index_l2x_Flrl_rofsub,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofi,i)   = (l2x_r%rAttr(index_l2x_Flrl_rofi,i) - &
                                               o2x_r%rAttr(index_o2x_Floo_qice,i))* lfrac !Jer: subtract off subshelf melt flux
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
            fldlist=seq_flds_x2r_fluxes, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
    end do
    call t_drvstopf  (trim(timer))

  end subroutine prep_rof_calc_l2r_rx

!================================================================================================

  subroutine prep_rof_calc_o2r_rx(fractions_ox, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2r_rx (note that o2r_rx is a local module variable)
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_ox(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri, eoi, efi
    character(*), parameter :: subname = '(prep_rof_calc_o2r_rx)'
    
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eoi = mod((eri-1),num_inst_ocn) + 1
       efi = mod((eri-1),num_inst_frc) + 1
       call seq_map_map(mapper_Fo2r, o2racc_ox(eoi), o2r_rx(eri), &
            fldlist='Floo_qice', norm=.true.)
    end do
    call t_drvstopf  (trim(timer))
  
  end subroutine  prep_rof_calc_o2r_rx !================================================================================================

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

  function prep_rof_get_o2racc_ox()
    type(mct_aVect), pointer :: prep_rof_get_o2racc_ox(:)
    prep_rof_get_o2racc_ox => o2racc_ox(:)
  end function prep_rof_get_o2racc_ox

  function prep_rof_get_o2racc_ox_cnt()
    integer, pointer :: prep_rof_get_o2racc_ox_cnt
    prep_rof_get_o2racc_ox_cnt => o2racc_ox_cnt
  end function prep_rof_get_o2racc_ox_cnt

  function prep_rof_get_mapper_Fo2r()
    type(seq_map), pointer :: prep_rof_get_mapper_Fo2r
    prep_rof_get_mapper_Fo2r => mapper_Fo2r
  end function prep_rof_get_mapper_Fo2r

end module prep_rof_mod
