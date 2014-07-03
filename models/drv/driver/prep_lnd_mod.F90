module prep_lnd_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8 
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_atm, num_inst_rof, num_inst_glc
  use seq_comm_mct    , only: num_inst_lnd, num_inst_frc
  use seq_comm_mct    , only: CPLID, LNDID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs                               
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata  
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: lnd, atm, rof, glc

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

  public :: prep_lnd_get_a2x_lx
  public :: prep_lnd_get_r2x_lx
  public :: prep_lnd_get_g2x_lx

  public :: prep_lnd_get_mapper_Sa2l
  public :: prep_lnd_get_mapper_Fa2l
  public :: prep_lnd_get_mapper_Fr2l
  public :: prep_lnd_get_mapper_SFg2l

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_lnd_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sa2l           ! needed in ccsm_comp_mod.F90 (setting of aream)
  type(seq_map), pointer :: mapper_Fa2l           ! needed in ccsm_comp_mod.F90 (seq_domain_check)
  type(seq_map), pointer :: mapper_Fr2l           ! needed in seq_frac_mct.F90
  type(seq_map), pointer :: mapper_SFg2l

  ! attribute vectors 
  type(mct_aVect), pointer :: a2x_lx(:) ! Atm export, lnd grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: r2x_lx(:) ! Rof export, lnd grid, lnd pes - allocated in lnd gc
  type(mct_aVect), pointer :: g2x_lx(:) ! Glc export, lnd grid, cpl pes - allocated in driver

  ! seq_comm_getData variables
  integer :: mpicom_CPLID                         ! MPI cpl communicator
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_lnd_init(infodata, atm_c2_lnd, rof_c2_lnd, glc_c2_lnd)
       
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
    !
    ! Local Variables
    integer                  :: lsize_l
    integer                  :: eai, eri, egi, eli
    logical                  :: samegrid_al   ! samegrid atm and land
    logical                  :: samegrid_lr   ! samegrid land and rof
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
    allocate(mapper_SFg2l)

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
          call mct_aVect_init(g2x_lx(egi), rList=seq_flds_g2x_fields, lsize=lsize_l)
          call mct_aVect_zero(g2x_lx(egi))
       end do

       samegrid_al = .true. 
       samegrid_lr = .true. 
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.

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
             write(logunit,F00) 'Initializing mapper_SFg2l'
          end if
          call seq_map_init_rearrolap(mapper_SFg2l, glc(1), lnd(1), 'mapper_SFg2l')
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_lnd_init

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
    integer                  :: eai, eri, egi, eli, efi
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
       call seq_map_map(mapper_SFg2l, g2x_gx, g2x_lx(egi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_g2x_lx

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

  function prep_lnd_get_mapper_SFg2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_SFg2l
    prep_lnd_get_mapper_SFg2l => mapper_SFg2l  
  end function prep_lnd_get_mapper_SFg2l

end module prep_lnd_mod
