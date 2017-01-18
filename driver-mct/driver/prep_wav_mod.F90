module prep_wav_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_atm, num_inst_ice, num_inst_ocn
  use seq_comm_mct    , only: num_inst_wav, num_inst_frc
  use seq_comm_mct    , only: CPLID, WAVID, logunit
  use seq_comm_mct    , only: seq_comm_getdata=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_getdata, seq_infodata_type
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: wav, ocn, ice, atm

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_wav_init
  public :: prep_wav_mrg

  public :: prep_wav_calc_a2x_wx
  public :: prep_wav_calc_o2x_wx
  public :: prep_wav_calc_i2x_wx

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_wav_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_sa2w
  type(seq_map), pointer :: mapper_so2w
  type(seq_map), pointer :: mapper_si2w

  ! attribute vectors
  type(mct_aVect), pointer :: o2x_wx(:) ! Ocn export, wav grid, cpl pes
  type(mct_aVect), pointer :: i2x_wx(:) ! Ice export, wav grid, cpl pes
  type(mct_aVect), pointer :: a2x_wx(:) ! Atm export, wav grid, cpl pes

  ! accumulation variables
  ! none at this time

  ! seq_comm_getData variables
  integer :: mpicom_CPLID                     ! MPI cpl communicator
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_wav_init(infodata, atm_c2_wav, ocn_c2_wav, ice_c2_wav)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_wav ! .true.  => atm to wav coupling on
    logical                 , intent(in)    :: ocn_c2_wav ! .true.  => ocn to wav coupling on
    logical                 , intent(in)    :: ice_c2_wav ! .true.  => ocn to wav coupling on
    !
    ! Local Variables
    integer                     :: eai , eoi, eii, ewi
    integer                     :: lsize_w
    logical                     :: samegrid_ow   ! samegrid ocean and wave
    logical                     :: samegrid_aw   ! samegrid atm and wave
    logical                     :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                     :: esmf_map_flag ! .true. => use esmf for mapping
    logical                     :: wav_present   ! .true. => wav is present
    character(CL)               :: atm_gnam      ! atm grid
    character(CL)               :: ocn_gnam      ! ocn grid
    character(CL)               :: wav_gnam      ! wav grid
    type(mct_avect) , pointer   :: w2x_wx
    character(*)    , parameter :: subname = '(prep_wav_init)'
    character(*)    , parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         wav_present=wav_present      , &
         ocn_gnam=ocn_gnam            , &
         wav_gnam=wav_gnam            , &
         atm_gnam=atm_gnam            , &
         esmf_map_flag=esmf_map_flag  )

    allocate(mapper_sa2w)
    allocate(mapper_so2w)
    allocate(mapper_si2w)

    if (wav_present) then

       call seq_comm_getData(CPLID, mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       w2x_wx => component_get_c2x_cx(wav(1))
       lsize_w = mct_aVect_lsize(w2x_wx)

       allocate(a2x_wx(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_wx(eai), rList=seq_flds_a2x_fields, lsize=lsize_w)
          call mct_aVect_zero(a2x_wx(eai))
       enddo
       allocate(o2x_wx(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_wx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_w)
          call mct_aVect_zero(o2x_wx(eoi))
       enddo
       allocate(i2x_wx(num_inst_ice))
       do eii = 1,num_inst_ice
          call mct_aVect_init(i2x_wx(eii), rList=seq_flds_i2x_fields, lsize=lsize_w)
          call mct_aVect_zero(i2x_wx(eii))
       enddo

       samegrid_ow = .true.
       samegrid_aw = .true.
       if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = .false.
       if (trim(atm_gnam) /= trim(wav_gnam)) samegrid_aw = .false.

       if (atm_c2_wav) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2w'
          end if
          call seq_map_init_rcfile(mapper_Sa2w, atm(1), wav(1), &
               'seq_maps.rc','atm2wav_smapname:','atm2wav_smaptype:',samegrid_aw, &
               'mapper_Sa2w initialization')
       endif
       call shr_sys_flush(logunit)
       if (ocn_c2_wav) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_So2w'
          end if
          call seq_map_init_rcfile(mapper_So2w, ocn(1), wav(1), &
               'seq_maps.rc','ocn2wav_smapname:','ocn2wav_smaptype:',samegrid_ow, &
               'mapper_So2w initialization')
       endif
       call shr_sys_flush(logunit)  !TODO ??? is this in Tony's code
       if (ice_c2_wav) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Si2w'
          end if
          call seq_map_init_rcfile(mapper_Si2w, ice(1), wav(1), &
               'seq_maps.rc','ice2wav_smapname:','ice2wav_smaptype:',samegrid_ow, &
               'mapper_Si2w initialization')
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_wav_init

  !================================================================================================

  subroutine prep_wav_mrg(infodata, fractions_wx, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge all wav inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_wx(:)
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eai, eoi, eii, ewi, efi
    type(mct_avect), pointer :: x2w_wx
    character(*), parameter  :: subname = '(prep_wav_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do ewi = 1,num_inst_wav
       ! Use fortran mod to address ensembles in merge
       eai = mod((ewi-1),num_inst_atm) + 1
       eoi = mod((ewi-1),num_inst_ocn) + 1
       eii = mod((ewi-1),num_inst_ice) + 1
       efi = mod((ewi-1),num_inst_frc) + 1

       x2w_wx => component_get_x2c_cx(wav(ewi))

       call prep_wav_merge(a2x_wx(eai), o2x_wx(eoi), i2x_wx(eii), fractions_wx(efi), x2w_wx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_wav_mrg

  !================================================================================================

  subroutine prep_wav_merge(a2x_w, o2x_w, i2x_w, frac_w, x2w_w)

    !-----------------------------------------------------------------------
    ! Arguments
    type(mct_aVect), intent(in)    :: a2x_w  ! input
    type(mct_aVect), intent(in)    :: o2x_w  ! input
    type(mct_aVect), intent(in)    :: i2x_w  ! input
    type(mct_aVect), intent(in)    :: frac_w ! input
    type(mct_aVect), intent(inout) :: x2w_w  ! output
    !-----------------------------------------------------------------------
    integer       :: nflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: o2x_sharedindices
    type(mct_aVect_sharedindices),save :: i2x_sharedindices
    character(*), parameter   :: subname = '(prep_wav_merge) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then
       nflds = mct_aVect_nRattr(x2w_w)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2w_w)
          mrgstr(i) = subname//'x2w%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(a2x_w, x2w_w, a2x_SharedIndices)
       ! QL, 150625, bug?
       ! a2x_SharedIndices -> o2x_SharedIndices
       ! a2x_SharedIndices -> i2x_SharedIndices
       call mct_aVect_setSharedIndices(o2x_w, x2w_w, o2x_SharedIndices)
       call mct_aVect_setSharedIndices(i2x_w, x2w_w, i2x_SharedIndices)

       !--- document copy operations ---
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, a2x_w)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
       enddo
       do i=1,o2x_SharedIndices%shared_real%num_indices
          i1=o2x_SharedIndices%shared_real%aVindices1(i)
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, o2x_w)
          mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field)
       enddo
       do i=1,i2x_SharedIndices%shared_real%num_indices
          i1=i2x_SharedIndices%shared_real%aVindices1(i)
          o1=i2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, i2x_w)
          mrgstr(o1) = trim(mrgstr(o1))//' = i2x%'//trim(field)
       enddo
    endif

    ! Create input wave state directly from atm, ocn, ice output state

    call mct_avect_zero(x2w_w)
    call mct_aVect_copy(aVin=a2x_w, aVout=x2w_w, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    call mct_aVect_copy(aVin=o2x_w, aVout=x2w_w, vector=mct_usevector, sharedIndices=o2x_SharedIndices)
    call mct_aVect_copy(aVin=i2x_w, aVout=x2w_w, vector=mct_usevector, sharedIndices=i2x_SharedIndices)

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

  end subroutine prep_wav_merge

  !================================================================================================

  subroutine prep_wav_calc_a2x_wx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create a2x_wx (note that a2x_wx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    type(mct_aVect), pointer :: a2x_ax
    character(*), parameter  :: subname = '(prep_wav_calc_a2x_wx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       a2x_ax => component_get_c2x_cx(atm(eai))
       call seq_map_map(mapper_Sa2w, a2x_ax, a2x_wx(eai), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_wav_calc_a2x_wx

  !================================================================================================

  subroutine prep_wav_calc_o2x_wx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_wx (note that o2x_wx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi
    type(mct_aVect), pointer :: o2x_ox
    character(*), parameter  :: subname = '(prep_wav_calc_o2x_wx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       o2x_ox => component_get_c2x_cx(ocn(eoi))
       call seq_map_map(mapper_So2w, o2x_ox, o2x_wx(eoi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_wav_calc_o2x_wx

  !================================================================================================

  subroutine prep_wav_calc_i2x_wx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create i2x_wx (note that i2x_wx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eii
    type(mct_aVect), pointer :: i2x_ix
    character(*), parameter   :: subname = '(prep_wav_calc_i2x_wx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eii = 1,num_inst_ice
       i2x_ix => component_get_c2x_cx(ice(eii))
       call seq_map_map(mapper_Si2w, i2x_ix, i2x_wx(eii), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_wav_calc_i2x_wx

end module prep_wav_mod
