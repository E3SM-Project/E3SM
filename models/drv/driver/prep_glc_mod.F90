module prep_glc_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8 
  use shr_kind_mod    , only: cs => SHR_KIND_CS
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
  public :: prep_glc_get_mapper_SFl2g

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_SFl2g

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

    allocate(mapper_SFl2g)

    if (glc_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)
       
       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)
       
       allocate(l2x_gx(num_inst_lnd))
       allocate(l2gacc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2x_gx(eli) ,lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))
          
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2gacc_lx(eli), lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       enddo
       l2gacc_lx_cnt = 0
       
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_SFl2g'
       end if
       call seq_map_init_rearrolap(mapper_SFl2g, lnd(1), glc(1), 'mapper_SFl2g')
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
  
  subroutine prep_glc_mrg(infodata, timer_mrg, timer_diag) 

    !---------------------------------------------------------------
    ! Description
    ! Merge glc inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    character(len=*)        , intent(in)    :: timer_mrg
    character(len=*)        , intent(in)    :: timer_diag
    !
    ! Local Variables
    integer :: egi, eli
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter  :: subname = '(prep_glc_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1

       x2g_gx => component_get_x2c_cx(glc(egi)) 
       call prep_glc_merge(l2x_gx(eli), x2g_gx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg

  !================================================================================================

  subroutine prep_glc_merge( s2x_g, x2g_g )

    !----------------------------------------------------------------------- 
    ! Arguments
    type(mct_aVect), intent(inout)  :: s2x_g  ! input
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !----------------------------------------------------------------------- 

    integer       :: nflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: s2x_sharedindices
    character(*), parameter   :: subname = '(prep_glc_merge) '

    !----------------------------------------------------------------------- 

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then
       nflds = mct_aVect_nRattr(x2g_g)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2g_g)
          mrgstr(i) = subname//'x2g%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(s2x_g, x2g_g, s2x_SharedIndices)

       !--- document copy operations ---
       do i=1,s2x_SharedIndices%shared_real%num_indices
          i1=s2x_SharedIndices%shared_real%aVindices1(i)
          o1=s2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, s2x_g)
          mrgstr(o1) = trim(mrgstr(o1))//' = s2x%'//trim(field)
       enddo
    endif

    ! Create input glc state directly from land snow output state
    call mct_aVect_copy(aVin=s2x_g, aVout=x2g_g, vector=mct_usevector, sharedIndices=s2x_SharedIndices)

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

  subroutine prep_glc_calc_l2x_gx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_gx (note that l2x_gx is a local module variable)
    ! Also l2x_gx is really the accumulated l2xacc_lx mapped to l2x_gx
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       call seq_map_map(mapper_SFl2g, l2gacc_lx(eli), l2x_gx(eli), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_l2x_gx

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

  function prep_glc_get_mapper_SFl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_SFl2g
    prep_glc_get_mapper_SFl2g => mapper_SFl2g  
  end function prep_glc_get_mapper_SFl2g

end module prep_glc_mod
