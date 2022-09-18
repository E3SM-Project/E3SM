module prep_lnd_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_kind_mod    , only: cxx => SHR_KIND_CXX
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_atm, num_inst_rof, num_inst_glc
  use seq_comm_mct    , only: num_inst_lnd, num_inst_frc
  use seq_comm_mct    , only: CPLID, LNDID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_comm_mct,     only: mlnid  ! iMOAB pid for land mesh on component pes
  use seq_comm_mct,     only: mhid     ! iMOAB id for atm instance
  use seq_comm_mct,     only: mphaid   ! iMOAB id for phys atm on atm pes
  use seq_comm_mct,     only: mhpgid   ! iMOAB id for atm pgx grid, on atm pes; created with se and gll grids
  use seq_comm_mct,     only: mblxid ! iMOAB id for mpas ocean migrated mesh to coupler pes
  use seq_comm_mct,     only: mbintxla ! iMOAB id for intx mesh between land and atmosphere
  use seq_comm_mct,     only: mbaxid   ! iMOAB id for atm migrated mesh to coupler pes
  use seq_comm_mct,     only: atm_pg_active  ! whether the atm uses FV mesh or not ; made true if fv_nphys > 0
  use dimensions_mod,   only: np     ! for atmosphere
  use seq_comm_mct,     only: seq_comm_getinfo => seq_comm_setptrs
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: lnd, atm, rof, glc
  use map_glc2lnd_mod   , only: map_glc2lnd_ec
  use iso_c_binding

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_lnd_init
  public :: prep_lnd_mrg

  public :: prep_lnd_calc_a2x_lx
  public :: prep_lnd_calc_r2x_lx
  public :: prep_lnd_calc_g2x_lx
  public :: prep_lnd_calc_z2x_lx

  public :: prep_lnd_get_a2x_lx
  public :: prep_lnd_get_r2x_lx
  public :: prep_lnd_get_g2x_lx
  public :: prep_lnd_get_z2x_lx

  public :: prep_lnd_get_mapper_Sa2l
  public :: prep_lnd_get_mapper_Fa2l
  public :: prep_lnd_get_mapper_Fr2l
  public :: prep_lnd_get_mapper_Sg2l
  public :: prep_lnd_get_mapper_Fg2l

  public :: prep_atm_lnd_moab ! it belongs here now

  public :: prep_lnd_migrate_moab
  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_lnd_merge
  private :: prep_lnd_set_glc2lnd_fields

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sa2l           ! needed in ccsm_comp_mod.F90 (setting of aream)
  type(seq_map), pointer :: mapper_Fa2l           ! needed in ccsm_comp_mod.F90 (seq_domain_check)
  type(seq_map), pointer :: mapper_Fr2l           ! needed in seq_frac_mct.F90
  type(seq_map), pointer :: mapper_Sg2l           ! currently unused (all g2l mappings use the flux mapper)
  type(seq_map), pointer :: mapper_Fg2l

  ! attribute vectors
  type(mct_aVect), pointer :: a2x_lx(:) ! Atm export, lnd grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: r2x_lx(:) ! Rof export, lnd grid, lnd pes - allocated in lnd gc
  type(mct_aVect), pointer :: g2x_lx(:) ! Glc export, lnd grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: z2x_lx(:) ! Iac export, lnd grid, cpl pes - allocated in driver

  ! seq_comm_getData variables
  integer :: mpicom_CPLID                         ! MPI cpl communicator

  ! field names and lists, for fields that need to be treated specially
  character(len=*), parameter :: glc_frac_field = 'Sg_ice_covered'
  character(len=*), parameter :: glc_topo_field = 'Sg_topo'
  character(len=*), parameter :: glc_icemask_field = 'Sg_icemask'
  ! fields mapped from glc to lnd, NOT separated by elevation class
  character(CXX) :: glc2lnd_non_ec_fields
  ! other fields (besides frac_field and topo_field) that are mapped from glc to lnd,
  ! separated by elevation class
  character(CXX) :: glc2lnd_ec_extra_fields
  !================================================================================================

#ifdef MOABDEBUG
  integer :: number_calls ! it is a static variable, used to count the number of projections
#endif
contains

  !================================================================================================

  subroutine prep_lnd_init(infodata, atm_c2_lnd, rof_c2_lnd, glc_c2_lnd, iac_c2_lnd)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_lnd ! .true.  => atm to lnd coupling on
    logical                 , intent(in)    :: rof_c2_lnd ! .true.  => rof to lnd coupling on
    logical                 , intent(in)    :: glc_c2_lnd ! .true.  => glc to lnd coupling on
    logical                 , intent(in)    :: iac_c2_lnd ! .true.  => iac to lnd coupling on
    !
    ! Local Variables
    integer                  :: lsize_l
    integer                  :: eai, eri, egi
    logical                  :: samegrid_al   ! samegrid atm and land
    logical                  :: samegrid_lr   ! samegrid land and rof
    logical                  :: samegrid_lg   ! samegrid land and glc
    logical                  :: esmf_map_flag ! .true. => use esmf for mapping
    logical                  :: lnd_present   ! .true. => land is present
    logical                  :: iamroot_CPLID ! .true. => CPLID masterproc
    character(CL)            :: atm_gnam      ! atm grid
    character(CL)            :: lnd_gnam      ! lnd grid
    character(CL)            :: rof_gnam      ! rof grid
    character(CL)            :: glc_gnam      ! glc grid
    type(mct_avect), pointer :: l2x_lx
    character(*), parameter  :: subname = '(prep_lnd_init)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         esmf_map_flag=esmf_map_flag,   &
         lnd_present=lnd_present,       &
         atm_gnam=atm_gnam,             &
         lnd_gnam=lnd_gnam,             &
         rof_gnam=rof_gnam,             &
         glc_gnam=glc_gnam)

    allocate(mapper_Sa2l)
    allocate(mapper_Fa2l)
    allocate(mapper_Fr2l)
    allocate(mapper_Sg2l)
    allocate(mapper_Fg2l)

    if (lnd_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       allocate(a2x_lx(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_lx(eai), rList=seq_flds_a2x_fields, lsize=lsize_l)
          call mct_aVect_zero(a2x_lx(eai))
       enddo
       allocate(r2x_lx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_aVect_init(r2x_lx(eri), rlist=seq_flds_r2x_fields, lsize=lsize_l)
          call mct_aVect_zero(r2x_lx(eri))
       end do
       allocate(g2x_lx(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(g2x_lx(egi), rList=seq_flds_x2l_fields_from_glc, lsize=lsize_l)
          call mct_aVect_zero(g2x_lx(egi))
       end do

       samegrid_al = .true.
       samegrid_lr = .true.
       samegrid_lg = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.
       if (trim(lnd_gnam) /= trim(glc_gnam)) samegrid_lg = .false.

       if (rof_c2_lnd) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fr2l'
          end if
          call seq_map_init_rcfile(mapper_Fr2l, rof(1), lnd(1), &
               'seq_maps.rc','rof2lnd_fmapname:','rof2lnd_fmaptype:',samegrid_lr, &
               string='mapper_Fr2l initialization',esmf_map=esmf_map_flag)
       end if
       call shr_sys_flush(logunit)

       if (atm_c2_lnd) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2l'
          end if
          call seq_map_init_rcfile(mapper_Sa2l, atm(1), lnd(1), &
               'seq_maps.rc','atm2lnd_smapname:','atm2lnd_smaptype:',samegrid_al, &
               'mapper_Sa2l initialization',esmf_map_flag)
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fa2l'
          end if
          call seq_map_init_rcfile(mapper_Fa2l, atm(1), lnd(1), &
               'seq_maps.rc','atm2lnd_fmapname:','atm2lnd_fmaptype:',samegrid_al, &
               'mapper_Fa2l initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)

       if (glc_c2_lnd) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sg2l'
          end if
          call seq_map_init_rcfile(mapper_Sg2l, glc(1), lnd(1), &
               'seq_maps.rc','glc2lnd_smapname:','glc2lnd_smaptype:',samegrid_lg, &
               'mapper_Sg2l initialization',esmf_map_flag)
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fg2l'
          end if
          call seq_map_init_rcfile(mapper_Fg2l, glc(1), lnd(1), &
               'seq_maps.rc','glc2lnd_fmapname:','glc2lnd_fmaptype:',samegrid_lg, &
               'mapper_Fg2l initialization',esmf_map_flag)

          call prep_lnd_set_glc2lnd_fields()
       endif
       call shr_sys_flush(logunit)

    end if
#ifdef MOABDEBUG
   number_calls = 0 ! it is a static variable, used to count the number of projections
#endif
  end subroutine prep_lnd_init

  !================================================================================================

  subroutine prep_lnd_set_glc2lnd_fields()

    !---------------------------------------------------------------
    ! Description
    ! Sets the module-level glc2lnd_non_ec_fields and glc2lnd_ec_extra_fields variables.
    !
    ! Local Variables
    character(len=CXX) :: temp_list

    character(*), parameter  :: subname = '(prep_lnd_set_glc2lnd_fields)'
    !---------------------------------------------------------------

    ! glc2lnd fields not separated by elevation class can be determined by finding fields
    ! that exist in both the g2x_to_lnd list and the x2l_from_glc list
    call shr_string_listIntersect(seq_flds_g2x_fields_to_lnd, &
         seq_flds_x2l_fields_from_glc, &
         glc2lnd_non_ec_fields)

    ! glc2lnd fields separated by elevation class are all fields not determined above.
    ! However, we also need to remove glc_frac_field and glc_topo_field from this list,
    ! because those are handled specially, so are not expected to be in this
    ! "extra_fields" list.
    !
    ! NOTE(wjs, 2015-04-24) I am going to the trouble of building this field list
    ! dynamically, rather than simply hard-coding the necessary fields (currently just
    ! 'Flgg_hflx'), so that new fields can be added in seq_flds_mod without needing to
    ! change any other code.
    call shr_string_listDiff(seq_flds_g2x_fields_to_lnd, &
         glc2lnd_non_ec_fields, &
         glc2lnd_ec_extra_fields)
    temp_list = glc2lnd_ec_extra_fields
    call shr_string_listDiff(temp_list, &
         glc_frac_field, &
         glc2lnd_ec_extra_fields)
    temp_list = glc2lnd_ec_extra_fields
    call shr_string_listDiff(temp_list, &
         glc_topo_field, &
         glc2lnd_ec_extra_fields)

  end subroutine prep_lnd_set_glc2lnd_fields

  !================================================================================================

  subroutine prep_lnd_mrg(infodata, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Prepare run phase, including running the merge
    !
    ! Arguments
    type(seq_infodata_type) , intent(in) :: infodata
    character(len=*)     , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eai, eri, egi, eli
    type(mct_aVect), pointer :: x2l_lx
    character(*), parameter  :: subname = '(prep_lnd_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       ! Use fortran mod to address ensembles in merge
       eai = mod((eli-1),num_inst_atm) + 1
       eri = mod((eli-1),num_inst_rof) + 1
       egi = mod((eli-1),num_inst_glc) + 1

       x2l_lx => component_get_x2c_cx(lnd(eli))  ! This is actually modifying x2l_lx
       call prep_lnd_merge( a2x_lx(eai), r2x_lx(eri), g2x_lx(egi), x2l_lx )
    enddo
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_lnd_mrg

  !================================================================================================

  subroutine prep_lnd_merge( a2x_l, r2x_l, g2x_l, x2l_l )
    !---------------------------------------------------------------
    ! Description
    ! Create input land state directly from atm, runoff and glc outputs
    !
    ! Arguments
    type(mct_aVect), intent(in)     :: a2x_l
    type(mct_aVect), intent(in)     :: r2x_l
    type(mct_aVect), intent(in)     :: g2x_l
    type(mct_aVect), intent(inout)  :: x2l_l
    !-----------------------------------------------------------------------
    integer       :: nflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: r2x_sharedindices
    type(mct_aVect_sharedindices),save :: g2x_sharedindices
    character(*), parameter   :: subname = '(prep_lnd_merge) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then
       nflds = mct_aVect_nRattr(x2l_l)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2l_l)
          mrgstr(i) = subname//'x2l%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(a2x_l, x2l_l, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(r2x_l, x2l_l, r2x_SharedIndices)
       call mct_aVect_setSharedIndices(g2x_l, x2l_l, g2x_SharedIndices)

       !--- document copy operations ---
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, a2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
       enddo
       do i=1,r2x_SharedIndices%shared_real%num_indices
          i1=r2x_SharedIndices%shared_real%aVindices1(i)
          o1=r2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, r2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = r2x%'//trim(field)
       enddo
       do i=1,g2x_SharedIndices%shared_real%num_indices
          i1=g2x_SharedIndices%shared_real%aVindices1(i)
          o1=g2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, g2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field)
       enddo
    endif

    call mct_aVect_copy(aVin=a2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    call mct_aVect_copy(aVin=r2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=r2x_SharedIndices)
    call mct_aVect_copy(aVin=g2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

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

  end subroutine prep_lnd_merge

  !================================================================================================

  subroutine prep_lnd_calc_a2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create  a2x_lx (note that a2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    type(mct_aVect), pointer :: a2x_ax
    character(*), parameter  :: subname = '(prep_lnd_calc_a2x_lx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       a2x_ax => component_get_c2x_cx(atm(eai))
       call seq_map_map(mapper_Fa2l, a2x_ax, a2x_lx(eai), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_a2x_lx

  !================================================================================================

  subroutine prep_lnd_calc_r2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create r2x_lx (note that r2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri
    type(mct_aVect) , pointer :: r2x_rx
    character(*), parameter :: subname = '(prep_lnd_calc_r2x_lx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       r2x_rx => component_get_c2x_cx(rof(eri))

       ! Note that one of these fields (a volr field) is remapped from rof -> lnd in
       ! map_lnd2rof_irrig_mod, because it is needed as a normalization term. So, if the
       ! details of this mapping call are changed in the future, it's possible that the
       ! equivalent r2l mapping in map_lnd2rof_irrig_mod should be changed to keep the two
       ! equivalent.
       call seq_map_map(mapper_Fr2l, r2x_rx, r2x_lx(eri), &
            fldlist=seq_flds_r2x_fluxes, norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_r2x_lx

  !================================================================================================

  subroutine prep_lnd_calc_g2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_lx (note that g2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_aVect), pointer :: g2x_gx
    character(*), parameter :: subname = '(prep_lnd_calc_g2x_lx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       g2x_gx => component_get_c2x_cx(glc(egi))

       ! Map fields that are NOT separated by elevation class on the land grid
       !
       ! These are mapped using a simple area-conservative remapping. (Note that we use
       ! the flux mapper even though these contain states, because we need these icemask
       ! fields to be mapped conservatively.)
       !
       ! Note that this mapping is redone for Sg_icemask in prep_glc_mod:
       ! prep_glc_map_qice_conservative_lnd2glc. If we ever change this mapping (e.g.,
       ! changing norm to .false.), then we should change the mapping there, too.
       !
       ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but
       ! this requires some more thought
       call seq_map_map(mapper_Fg2l, g2x_gx, g2x_lx(egi), &
            fldlist = glc2lnd_non_ec_fields, norm=.true.)

       ! Map fields that are separated by elevation class on the land grid
       call map_glc2lnd_ec( &
            g2x_g = g2x_gx, &
            frac_field = glc_frac_field, &
            topo_field = glc_topo_field, &
            icemask_field = glc_icemask_field, &
            extra_fields = glc2lnd_ec_extra_fields, &
            mapper = mapper_Fg2l, &
            g2x_l = g2x_lx(egi))
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_g2x_lx

  !================================================================================================

  subroutine prep_lnd_calc_z2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create z2x_lx (note that z2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_aVect), pointer :: z2x_gx
    character(*), parameter :: subname = '(prep_lnd_calc_z2x_lx)'
    !---------------------------------------------------------------

    ! Stub

  end subroutine prep_lnd_calc_z2x_lx

  !================================================================================================

  function prep_lnd_get_a2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_a2x_lx(:)
    prep_lnd_get_a2x_lx => a2x_lx(:)
  end function prep_lnd_get_a2x_lx

  function prep_lnd_get_r2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_r2x_lx(:)
    prep_lnd_get_r2x_lx => r2x_lx(:)
  end function prep_lnd_get_r2x_lx

  function prep_lnd_get_g2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_g2x_lx(:)
    prep_lnd_get_g2x_lx => g2x_lx(:)
  end function prep_lnd_get_g2x_lx

  function prep_lnd_get_z2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_z2x_lx(:)
    prep_lnd_get_z2x_lx => z2x_lx(:)
  end function prep_lnd_get_z2x_lx

  function prep_lnd_get_mapper_Sa2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Sa2l
    prep_lnd_get_mapper_Sa2l => mapper_Sa2l
  end function prep_lnd_get_mapper_Sa2l

  function prep_lnd_get_mapper_Fa2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Fa2l
    prep_lnd_get_mapper_Fa2l => mapper_Fa2l
  end function prep_lnd_get_mapper_Fa2l

  function prep_lnd_get_mapper_Fr2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Fr2l
    prep_lnd_get_mapper_Fr2l => mapper_Fr2l
  end function prep_lnd_get_mapper_Fr2l

  function prep_lnd_get_mapper_Sg2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Sg2l
    prep_lnd_get_mapper_Sg2l => mapper_Sg2l
  end function prep_lnd_get_mapper_Sg2l

  function prep_lnd_get_mapper_Fg2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Fg2l
    prep_lnd_get_mapper_Fg2l => mapper_Fg2l
  end function prep_lnd_get_mapper_Fg2l

  ! moved from prep_atm
  subroutine prep_atm_lnd_moab(infodata)

   use iMOAB, only: iMOAB_CoverageGraph, iMOAB_ComputeScalarProjectionWeights, iMOAB_ComputeCommGraph
   use iMOAB, only: iMOAB_DefineTagStorage
   !---------------------------------------------------------------
   ! Description
   ! If the land is on the same mesh as atm, we do not need to compute intx
   !  Just use compute graph between phys atm and lnd on coupler, to be able to send
   !  data from atm phys to atm on coupler for projection on land
   ! in the tri-grid case, atm and land use different meshes, so use coverage anyway
   !
   ! Arguments
   type(seq_infodata_type) , intent(in)    :: infodata

   character(*), parameter          :: subname = '(prep_atm_lnd_moab)'
   integer :: ierr

   logical                          :: atm_present    ! .true.  => atm is present
   logical                          :: lnd_present    ! .true.  => lnd is present
   integer                  :: id_join
   integer                  :: mpicom_join
   integer                  :: context_id ! used to define context for coverage (this case, land on coupler)
   integer                  :: atm_id
   character*32             :: dm1, dm2, dofnameATM, dofnameLND, wgtIdef
   integer                  :: orderLND, orderATM, volumetric, fInverseDistanceMap, noConserve, validate
   integer                  :: fNoBubble, monotonicity
   integer                  :: mpigrp_CPLID ! coupler pes group, used for comm graph phys <-> atm-ocn
   integer                  :: mpigrp_old   !  component group pes (phys grid atm) == atm group
   integer                  :: typeA, typeB ! type for computing graph;
   integer                  :: idintx ! in this case, id of moab intersection between atm and lnd, on coupler pes
                                      ! used only for tri-grid case
   integer                  :: tagtype, numco, tagindex
   character(CXX)            :: tagname    ! will store all seq_flds_a2x_fields
   character(CL)            :: atm_gnam      ! atm grid
   character(CL)            :: lnd_gnam      ! lnd grid
   logical                  :: samegrid_al 

   call seq_infodata_getData(infodata, &
        atm_present=atm_present,       &
        lnd_present=lnd_present,       &
         atm_gnam=atm_gnam,             &
         lnd_gnam=lnd_gnam)

    samegrid_al = .true.
    if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
   !  it involves initial atm app; mhid; or pg2 mesh , in case atm_pg_active also migrate atm mesh on coupler pes, mbaxid
   !  intx lnd atm are in mbintxla ; remapper also has some info about coverage mesh
   ! after this, the sending of tags from atm pes to coupler pes, in land context will use the new par
   ! comm graph, that has more precise info about
   ! how to get mpicomm for joint atm + coupler
   id_join = atm(1)%cplcompid
   atm_id   = atm(1)%compid
   ! maybe we can use a moab-only id, defined like mbintxao, mhid, somewhere else (seq_comm_mct)
   ! we cannot use mbintxla because it may not exist on atm comp yet;
   context_id = lnd(1)%cplcompid
   call seq_comm_getinfo(ID_join,mpicom=mpicom_join)
   if ( .not. samegrid_al ) then
     if (atm_pg_active ) then ! use mhpgid mesh
       ierr = iMOAB_CoverageGraph(mpicom_join, mhpgid, mbaxid, mbintxla, atm_id, id_join, context_id);
     else
       ierr = iMOAB_CoverageGraph(mpicom_join, mhid, mbaxid, mbintxla, atm_id, id_join, context_id);
     endif
   else
     ! this is the moment we compute the comm graph between phys grid atm and land on coupler pes.
     ! We do not need to compute intersection in this case, as the DOFs are exactly the same
     !  see imoab_phatm_ocn_coupler.cpp in MOAB source code, no intx needed, just compute graph
     typeA = 2 ! point cloud 
     typeB = 2 ! 
     call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)   !  second group, the coupler group CPLID is global variable
     call seq_comm_getinfo(atm_id, mpigrp=mpigrp_old) 
     !  context_id = lnd(1)%cplcompid
     ierr = iMOAB_ComputeCommGraph( mphaid, mblxid, mpicom_join, mpigrp_old, mpigrp_CPLID, &
       typeA, typeB, atm_id, context_id)
     if (ierr .ne. 0) then
       write(logunit,*) subname,' error in computing graph phys grid - lnd on coupler '
       call shr_sys_abort(subname//' ERROR  in computing graph phys grid - lnd on coupler ')
     endif

   endif
   if (ierr .ne. 0) then
     write(logunit,*) subname,' error in computing coverage graph atm/lnd '
     call shr_sys_abort(subname//' ERROR in computing coverage graph atm/lnd ')
   endif

   ! this is true only for tri-grid cases
   if (mbintxla .ge. 0 ) then ! weights are computed over coupler pes
     ! copy from atm - ocn  , it is now similar, as land is full mesh, not pc cloud
     wgtIdef = 'scalar'//C_NULL_CHAR
     volumetric = 0 !  TODO: check this , for PC ; for imoab_coupler test, volumetric is 0
     if (atm_pg_active) then
       dm1 = "fv"//C_NULL_CHAR
       dofnameATM="GLOBAL_ID"//C_NULL_CHAR
       orderATM = 1 !  fv-fv
     else
       dm1 = "cgll"//C_NULL_CHAR
       dofnameATM="GLOBAL_DOFS"//C_NULL_CHAR
       orderATM = np !  it should be 4
       volumetric = 1 
     endif

     dofnameLND="GLOBAL_ID"//C_NULL_CHAR
     orderLND = 1  !  not much arguing
     
     ! is the land mesh explicit or point cloud ? based on samegrid_al flag:
     if (samegrid_al) then
       dm2 = "pcloud"//C_NULL_CHAR
       wgtIdef = 'scalar-pc'//C_NULL_CHAR
     else
       dm2 = "fv"//C_NULL_CHAR ! land is FV
     endif
     fNoBubble = 1
     monotonicity = 0 !
     noConserve = 0
     validate = 0
     fInverseDistanceMap = 0

     ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxla, wgtIdef, &
                                               trim(dm1), orderATM, trim(dm2), orderLND, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameATM), trim(dofnameLND) )
     if (ierr .ne. 0) then
       write(logunit,*) subname,' error in computing weights atm land '
       call shr_sys_abort(subname//' ERROR  in computing weights atm land')
     endif
   endif
   ! we will use intx atm-lnd mesh only when land is explicit
   if (.not. samegrid_al) then
     ! as with ocn, data is sent from atm ph to the intx atm/lnd, not from pg2 mesh anymore
     ! for that, we will use the comm graph between atm ph and atm pg2 intersected with land!
     ! copy from ocn logic, just replace with land
         ! compute the comm graph between phys atm and intx-atm-lnd, to be able to send directly from phys atm
     ! towards coverage mesh on atm for intx to land / now that land is full mesh!
     ! this is similar to imoab_phatm_ocn_coupler.cpp test in moab
     !    int typeA = 2; // point cloud
     call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)   !  second group, the coupler group CPLID is global variable
     call seq_comm_getinfo(atm_id, mpigrp=mpigrp_old)    !  component group pes, from atm id ( also ATMID(1) )

     typeA = 2 ! point cloud, phys atm in this case
     ! idintx is a unique number of MOAB app that takes care of intx between lnd and atm mesh
     idintx = 100*atm(1)%cplcompid + lnd(1)%cplcompid ! something different, to differentiate it; ~ 600+lnd !
     if (atm_pg_active) then
       typeB = 3 ! fv on atm side too !! imoab_apg2_ol  coupler example
                 ! atm cells involved in intersection (pg 2 in this case)
                 ! this will be used now to send
                 ! data from phys grid directly to atm-lnd intx , for later projection
                 ! context is the same, atm - lnd intx id !

     else
       typeB = 1 ! atm cells involved in intersection (spectral in this case) ! this will be used now to send
                 ! data from phys grid directly to atm-lnd intx , for later projection
                 ! context is the same, atm - lnd intx id !
     endif
     ierr = iMOAB_ComputeCommGraph( mphaid, mbintxla, mpicom_join, mpigrp_old, mpigrp_CPLID, &
           typeA, typeB, atm_id, idintx)
     if (ierr .ne. 0) then
       write(logunit,*) subname,' error in computing graph phys grid - atm/lnd intx '
       call shr_sys_abort(subname//' ERROR  in computing graph phys grid - atm/ocn intx ')
     endif
   endif ! if (.not. samegrid_al)

   if (mblxid .ge. 0) then
      ! in any case, we need to define the tags on landx from the phys atm seq_flds_a2x_fields
      tagtype = 1  ! dense, double
      numco = 1 !  one value per vertex / entity
      tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
      ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
      if ( ierr > 0) then
         call shr_sys_abort(subname//' fail to define seq_flds_a2x_fields for lnd x moab mesh ')
      endif
   endif

 end subroutine prep_atm_lnd_moab

  ! exposed method to migrate projected tag from coupler pes back to land pes
  subroutine prep_lnd_migrate_moab(infodata)

   use iMOAB, only: iMOAB_SendElementTag, iMOAB_ReceiveElementTag, iMOAB_FreeSenderBuffers, &
      iMOAB_WriteMesh
  !---------------------------------------------------------------
    ! Description
    ! After a2lTbot_proj, a2lVbot_proj, a2lUbot_proj were computed on lnd mesh on coupler, they need
    !   to be migrated to the land pes
    !  maybe the land solver will use it (later)?
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata

    integer :: ierr

    logical                          :: atm_present    ! .true.  => atm is present
    logical                          :: lnd_present    ! .true.  => lnd is present
    integer                  :: id_join
    integer                  :: mpicom_join
    integer                  :: lndid1
    integer                  :: context_id
    character*32             :: dm1, dm2
    character*50             :: tagName
    character*32             :: outfile, wopts, lnum
    integer                  :: orderLND, orderATM, volumetric, noConserve, validate 

    call seq_infodata_getData(infodata, &
         atm_present=atm_present,       &
         lnd_present=lnd_present)

  !  it involves initial ocn app; mpoid; also migrated ocn mesh mesh on coupler pes, mbaxid
  ! after this, the sending of tags from coupler pes to ocn pes will use initial graph
       !  (not processed for coverage)
  ! how to get mpicomm for joint ocn + coupler
    id_join = lnd(1)%cplcompid
    lndid1   = lnd(1)%compid
!     call seq_comm_getinfo(ID_join,mpicom=mpicom_join)
!     context_id = -1
!     ! now send the tag a2oTbot_proj, a2oUbot_proj, a2oVbot_proj from ocn on coupler pes towards original ocean mesh
!     tagName = 'a2lTbot_proj:a2lUbot_proj:a2lVbot_proj:'//C_NULL_CHAR !  defined in prep_atm_mod.F90!!!

!     if (mblxid .ge. 0) then !  send because we are on coupler pes

!       ! basically, use the initial partitioning
!        context_id = lndid1
!        ierr = iMOAB_SendElementTag(mblxid, tagName, mpicom_join, context_id)

!     endif
!     if (mlnid .ge. 0 ) then !  we are on land pes, for sure
!       ! receive on land pes, a tag that was computed on coupler pes
!        context_id = id_join
!        ierr = iMOAB_ReceiveElementTag(mlnid, tagName, mpicom_join, context_id)
!     !CHECKRC(ierr, "cannot receive tag values")
!     endif

!     ! we can now free the sender buffers
!     if (mblxid .ge. 0) then
!        context_id = lndid1
!        ierr = iMOAB_FreeSenderBuffers(mblxid, context_id)
!        ! CHECKRC(ierr, "cannot free buffers used to send projected tag towards the ocean mesh")
!     endif

! #ifdef MOABDEBUG
!     if (mlnid .ge. 0 ) then !  we are on land pes, for sure
!       number_calls = number_calls + 1
!       write(lnum,"(I0.2)") number_calls
!       outfile = 'wholeLND_proj'//trim(lnum)//'.h5m'//C_NULL_CHAR
!       wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
!       ierr = iMOAB_WriteMesh(mlnid, trim(outfile), trim(wopts))
!     endif
! #endif

  end subroutine prep_lnd_migrate_moab

end module prep_lnd_mod
