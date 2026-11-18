module prep_wav_mod
   use cplcomp_moab_helpers_mod, only: moab_mesh_discretization

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cs => SHR_KIND_CS
   use shr_kind_mod    , only: cl => SHR_KIND_CL
   use shr_kind_mod    , only: CXX => SHR_KIND_CXX
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use shr_string_mod  , only: shr_string_listIntersect, shr_string_listGetIndexF
  use seq_comm_mct    , only: num_inst_atm, num_inst_ice, num_inst_ocn
  use seq_comm_mct    , only: num_inst_wav, num_inst_frc
  use seq_comm_mct    , only: CPLID, WAVID, logunit
  use seq_comm_mct    , only: seq_comm_getdata=>seq_comm_setptrs
  use seq_comm_mct    , only: mbaxid, mboxid, mbixid
  use seq_comm_mct    , only: mbwxid
  use seq_comm_mct    , only: mbintxaw, mbintxow, mbintxiw
  use seq_infodata_mod, only: seq_infodata_getdata, seq_infodata_type
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: wav, ocn, ice, atm
  use iso_c_binding

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_wav_init
  public :: prep_wav_mrg_moab

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

  ! atm/ocn/ice exports WW3 imports, as in the MCT x2w merge; set in prep_wav_init
  character(CXX) :: a2w_states_fldlist = ''
  character(CXX) :: o2w_states_fldlist = ''
  character(CXX) :: i2w_states_fldlist = ''

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

      use iMOAB, only: iMOAB_RegisterApplication, iMOAB_ComputeCommGraph

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_wav ! .true.  => atm to wav coupling on
    logical                 , intent(in)    :: ocn_c2_wav ! .true.  => ocn to wav coupling on
    logical                 , intent(in)    :: ice_c2_wav ! .true.  => ice to wav coupling on
    !
    ! Local Variables
    integer                     :: eai , eoi, eii
    integer                     :: lsize_w
    logical                     :: samegrid_ow   ! samegrid ocean and wave
    logical                     :: samegrid_aw   ! samegrid atm and wave
    logical                     :: samegrid_iw   ! samegrid ice and wave
    logical                     :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                     :: esmf_map_flag ! .true. => use esmf for mapping
    logical                     :: wav_present   ! .true. => wav is present
    character(CL)               :: atm_gnam      ! atm grid
    character(CL)               :: ocn_gnam      ! ocn grid
    character(CL)               :: ice_gnam      ! ice grid
    character(CL)               :: wav_gnam      ! wav grid
    type(mct_avect) , pointer   :: w2x_wx
   integer                     :: ierr, idintx
   integer                     :: type1, arearead, src_type, mpigrp_CPLID
    character*32                :: appname
    character*32                :: wgtIdSa2w, wgtIdSo2w, wgtIdSi2w
    character(*)    , parameter :: subname = '(prep_wav_init)'
    character(*)    , parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         wav_present=wav_present      , &
         ocn_gnam=ocn_gnam            , &
         ice_gnam=ice_gnam            , &
         wav_gnam=wav_gnam            , &
         atm_gnam=atm_gnam            , &
         esmf_map_flag=esmf_map_flag  )

    wgtIdSa2w = 'scalar_a2w'
    wgtIdSo2w = 'scalar_o2w'
    wgtIdSi2w = 'scalar_i2w'

    allocate(mapper_sa2w)
    allocate(mapper_so2w)
    allocate(mapper_si2w)

    if (wav_present) then

       call seq_comm_getData(CPLID, mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)
      call seq_comm_getdata(CPLID, mpigrp=mpigrp_CPLID)

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
       samegrid_iw = .true.
       if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = .false.
       if (trim(atm_gnam) /= trim(wav_gnam)) samegrid_aw = .false.
       if (trim(ice_gnam) /= trim(wav_gnam)) samegrid_iw = .false.

      if (atm_c2_wav .and. wav_atm_coup /= 'none') then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2w'
          end if
          call seq_map_mapinit(mapper_Sa2w, mpicom_CPLID)
          call shr_string_listIntersect(seq_flds_a2x_fields, seq_flds_x2w_fields, a2w_states_fldlist)
          if ( (mbaxid .ge. 0) .and. (mbwxid .ge. 0) ) then
             mapper_Sa2w%src_mbid = mbaxid
             mapper_Sa2w%tgt_mbid = mbwxid
             mapper_Sa2w%src_context = atm(1)%cplcompid
             mapper_Sa2w%weight_identifier = wgtIdSa2w
             mapper_Sa2w%mbname = 'mapper_Sa2w'
             if (samegrid_aw) then
                mapper_Sa2w%rearrange_only = .true.
                mapper_Sa2w%strategy = "rearrange"
                mapper_Sa2w%intx_context = wav(1)%cplcompid
                src_type = moab_mesh_discretization(mbaxid)
                ierr = iMOAB_ComputeCommGraph(mbaxid, mbwxid, mpicom_CPLID, &
                     mpigrp_CPLID, mpigrp_CPLID, src_type, 2, atm(1)%cplcompid, wav(1)%cplcompid)
                if (ierr /= 0) call shr_sys_abort(subname//' cannot create same-grid ATM-WAV graph')
             else
                if (mbintxaw < 0) then
                   appname = "ATM_WAV_COU"
                   idintx = 100*atm(1)%cplcompid + wav(1)%cplcompid
                   ierr = iMOAB_RegisterApplication(trim(appname)//C_NULL_CHAR, mpicom_CPLID, idintx, mbintxaw)
                   if (ierr .ne. 0) then
                      write(logunit,*) subname,' error in registering atm-wave intersection'
                      call shr_sys_abort(subname//' ERROR in registering atm-wave intersection')
                   end if
                end if
                mapper_Sa2w%intx_mbid = mbintxaw
                mapper_Sa2w%intx_context = 100*atm(1)%cplcompid + wav(1)%cplcompid
                type1 = 3
                arearead = 0
                call moab_map_init_rcfile(mapper_Sa2w, type1, &
                     'seq_maps.rc', 'atm2wav_smapname:', 'atm2wav_smaptype:', samegrid_aw, &
                   arearead, wgtIdSa2w, 'mapper_Sa2w MOAB initialization', esmf_map_flag, &
                   source_discretization_type=moab_mesh_discretization(mbaxid), &
                   target_discretization_type=2)
                   call seq_map_setup_coverage(mapper_Sa2w)
             endif
          endif
       endif
       call shr_sys_flush(logunit)
      if (ocn_c2_wav .and. wav_ocn_coup /= 'none') then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_So2w'
          end if
          call seq_map_mapinit(mapper_So2w, mpicom_CPLID)
          call shr_string_listIntersect(seq_flds_o2x_fields, seq_flds_x2w_fields, o2w_states_fldlist)
          if ( (mboxid .ge. 0) .and. (mbwxid .ge. 0) ) then
             mapper_So2w%src_mbid = mboxid
             mapper_So2w%tgt_mbid = mbwxid
             mapper_So2w%src_context = ocn(1)%cplcompid
             mapper_So2w%weight_identifier = wgtIdSo2w
             mapper_So2w%mbname = 'mapper_So2w'
             if (samegrid_ow) then
                mapper_So2w%rearrange_only = .true.
                mapper_So2w%strategy = "rearrange"
                mapper_So2w%intx_context = wav(1)%cplcompid
                src_type = moab_mesh_discretization(mboxid)
                ierr = iMOAB_ComputeCommGraph(mboxid, mbwxid, mpicom_CPLID, &
                     mpigrp_CPLID, mpigrp_CPLID, src_type, 2, ocn(1)%cplcompid, wav(1)%cplcompid)
                if (ierr /= 0) call shr_sys_abort(subname//' cannot create same-grid OCN-WAV graph')
             else
                if (mbintxow < 0) then
                   appname = "OCN_WAV_COU"
                   idintx = 100*ocn(1)%cplcompid + wav(1)%cplcompid
                   ierr = iMOAB_RegisterApplication(trim(appname)//C_NULL_CHAR, mpicom_CPLID, idintx, mbintxow)
                   if (ierr .ne. 0) then
                      write(logunit,*) subname,' error in registering ocn-wave intersection'
                      call shr_sys_abort(subname//' ERROR in registering ocn-wave intersection')
                   end if
                end if
                mapper_So2w%intx_mbid = mbintxow
                mapper_So2w%intx_context = 100*ocn(1)%cplcompid + wav(1)%cplcompid
                type1 = 3
                arearead = 0
                call moab_map_init_rcfile(mapper_So2w, type1, &
                     'seq_maps.rc', 'ocn2wav_smapname:', 'ocn2wav_smaptype:', samegrid_ow, &
                   arearead, wgtIdSo2w, 'mapper_So2w MOAB initialization', esmf_map_flag, &
                   source_discretization_type=moab_mesh_discretization(mboxid), &
                   target_discretization_type=2)
                   call seq_map_setup_coverage(mapper_So2w)
             endif
          endif
       endif
       call shr_sys_flush(logunit)  !TODO ??? is this in Tony's code
      if (ice_c2_wav .and. wav_ice_coup /= 'none') then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Si2w'
          end if
          call seq_map_mapinit(mapper_Si2w, mpicom_CPLID)
          ! Ice states WW3 imports; Si_ifloe is included only for two-way coupling
          call shr_string_listIntersect(seq_flds_i2x_fields, seq_flds_x2w_fields, i2w_states_fldlist)
          ! WW3 needs the ice fraction for both one-way and two-way coupling
          if (shr_string_listGetIndexF(i2w_states_fldlist, 'Si_ifrac') == 0) then
             call shr_sys_abort(subname//' ERROR Si_ifrac missing from ice to wave fields')
          endif
          if ( (mbixid .ge. 0) .and. (mbwxid .ge. 0) ) then
             mapper_Si2w%src_mbid = mbixid
             mapper_Si2w%tgt_mbid = mbwxid
             mapper_Si2w%src_context = ice(1)%cplcompid
             mapper_Si2w%weight_identifier = wgtIdSi2w
             mapper_Si2w%mbname = 'mapper_Si2w'
             if (samegrid_iw) then
                mapper_Si2w%rearrange_only = .true.
                mapper_Si2w%strategy = "rearrange"
                mapper_Si2w%intx_context = wav(1)%cplcompid
                src_type = moab_mesh_discretization(mbixid)
                ierr = iMOAB_ComputeCommGraph(mbixid, mbwxid, mpicom_CPLID, &
                     mpigrp_CPLID, mpigrp_CPLID, src_type, 2, ice(1)%cplcompid, wav(1)%cplcompid)
                if (ierr /= 0) call shr_sys_abort(subname//' cannot create same-grid ICE-WAV graph')
             else
                if (mbintxiw < 0) then
                   appname = "ICE_WAV_COU"
                   idintx = 100*ice(1)%cplcompid + wav(1)%cplcompid
                   ierr = iMOAB_RegisterApplication(trim(appname)//C_NULL_CHAR, mpicom_CPLID, idintx, mbintxiw)
                   if (ierr .ne. 0) then
                      write(logunit,*) subname,' error in registering ice-wave intersection'
                      call shr_sys_abort(subname//' ERROR in registering ice-wave intersection')
                   end if
                end if
                mapper_Si2w%intx_mbid = mbintxiw
                mapper_Si2w%intx_context = 100*ice(1)%cplcompid + wav(1)%cplcompid
                type1 = 3
                arearead = 0
                call moab_map_init_rcfile(mapper_Si2w, type1, &
                     'seq_maps.rc', 'ice2wav_smapname:', 'ice2wav_smaptype:', samegrid_iw, &
                   arearead, wgtIdSi2w, 'mapper_Si2w MOAB initialization', esmf_map_flag, &
                   source_discretization_type=moab_mesh_discretization(mbixid), &
                   target_discretization_type=2)
                   call seq_map_setup_coverage(mapper_Si2w)
             endif
          endif
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_wav_init

  !================================================================================================

  subroutine prep_wav_mrg_moab(infodata, fractions_wx, timer_mrg)

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
    character(*), parameter  :: subname = '(prep_wav_mrg_moab)'
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

   end subroutine prep_wav_mrg_moab

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
       if (len_trim(a2w_states_fldlist) > 0) &
            call seq_map_map(mapper_Sa2w, a2x_ax, a2x_wx(eai), fldlist=trim(a2w_states_fldlist), norm=.true.)
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
       if (len_trim(o2w_states_fldlist) > 0) &
            call seq_map_map(mapper_So2w, o2x_ox, o2x_wx(eoi), fldlist=trim(o2w_states_fldlist), norm=.true.)
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
       if (len_trim(i2w_states_fldlist) > 0) &
            call seq_map_map(mapper_Si2w, i2x_ix, i2x_wx(eii), fldlist=trim(i2w_states_fldlist), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_wav_calc_i2x_wx

end module prep_wav_mod
