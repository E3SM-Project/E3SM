module dice_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_DistGrid, ESMF_DistGridGet
  use ESMF                  , only : ESMF_State, ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF                  , only : operator(/=), operator(==)
  use perf_mod              , only : t_startf, t_stopf, t_barrierf
  use mct_mod               , only : mct_gsmap_init, mct_avect_indexRA 
  use shr_kind_mod          , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod         , only : shr_const_pi, shr_const_spval, shr_const_tkfrz, shr_const_latice
  use shr_mpi_mod           , only : shr_mpi_bcast
  use shr_frz_mod           , only : shr_frz_freezetemp
  use shr_cal_mod           , only : shr_cal_calendarname, shr_cal_datetod2string
  use shr_sys_mod           , only : shr_sys_abort
  use shr_strdata_mod       , only : shr_strdata_init_model_domain
  use shr_strdata_mod       , only : shr_strdata_init_streams
  use shr_strdata_mod       , only : shr_strdata_init_mapping
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_pioinit
  use shr_strdata_mod       , only : shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod       , only : shr_strdata_advance, shr_strdata_restWrite
  use dshr_methods_mod      , only : chkerr, state_getfldptr
  use dshr_nuopc_mod        , only : fld_list_type, dshr_fld_add
  use dice_flux_atmice_mod  , only : dice_flux_atmice
  use dice_shr_mod          , only : datamode       ! namelist input
  use dice_shr_mod          , only : rest_file      ! namelist input
  use dice_shr_mod          , only : rest_file_strm ! namelist input
  use dice_shr_mod          , only : flux_swpf      ! namelist input -short-wave penatration factor
  use dice_shr_mod          , only : flux_Qmin      ! namelist input -bound on melt rate
  use dice_shr_mod          , only : flux_Qacc      ! namelist input -activates water accumulation/melt wrt Q
  use dice_shr_mod          , only : flux_Qacc0     ! namelist input -initial water accumulation value
  use dice_shr_mod          , only : nullstr
  use dice_shr_mod          , only : SDICE
  use shr_pcdf_mod

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dice_comp_advertise
  public :: dice_comp_init
  public :: dice_comp_run
  public :: dice_comp_setptrs

  public :: dice_comp_copy_streams

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  integer                    :: debug_import = 0      ! debug level (if > 0 will print all import fields)
  integer                    :: debug_export = 0      ! debug level (if > 0 will print all export fields)

  real(R8),parameter         :: pi     = shr_const_pi      ! pi
  real(R8),parameter         :: spval  = shr_const_spval   ! flags invalid data
  real(R8),parameter         :: tFrz   = shr_const_tkfrz   ! temp of freezing
  real(R8),parameter         :: latice = shr_const_latice  ! latent heat of fusion
  real(R8),parameter         :: waterMax = 1000.0_R8       ! wrt iFrac comp & frazil ice (kg/m^2)

  !----- surface albedo constants ------
  real(R8),parameter         :: snwfrac = 0.286_R8 ! snow cover fraction ~ [0,1]
  real(R8),parameter         :: as_nidf = 0.950_R8 ! albedo: snow,near-infr,diffuse
  real(R8),parameter         :: as_vsdf = 0.700_R8 ! albedo: snow,visible  ,diffuse
  real(R8),parameter         :: as_nidr = 0.960_R8 ! albedo: snow,near-infr,direct
  real(R8),parameter         :: as_vsdr = 0.800_R8 ! albedo: snow,visible  ,direct
  real(R8),parameter         :: ai_nidf = 0.700_R8 ! albedo: ice, near-infr,diffuse
  real(R8),parameter         :: ai_vsdf = 0.500_R8 ! albedo: ice, visible  ,diffuse
  real(R8),parameter         :: ai_nidr = 0.700_R8 ! albedo: ice, near-infr,direct
  real(R8),parameter         :: ai_vsdr = 0.500_R8 ! albedo: ice, visible  ,direct
  real(R8),parameter         :: ax_nidf = ai_nidf*(1.0_R8-snwfrac) + as_nidf*snwfrac
  real(R8),parameter         :: ax_vsdf = ai_vsdf*(1.0_R8-snwfrac) + as_vsdf*snwfrac
  real(R8),parameter         :: ax_nidr = ai_nidr*(1.0_R8-snwfrac) + as_nidr*snwfrac
  real(R8),parameter         :: ax_vsdr = ai_vsdr*(1.0_R8-snwfrac) + as_vsdr*snwfrac

  integer                    :: index_lat, index_lon
  integer     , pointer      :: imask(:)
  real(R8)    , pointer      :: xc(:), yc(:)       ! arrays of model latitudes and longitudes
  real(R8)    , pointer      :: water(:)
  real(R8)    , pointer      :: tfreeze(:)
  !real(R8)   , pointer      :: ifrac0(:)

  logical                    :: firstcall = .true. ! first call logical
  character(len=*),parameter :: rpfile = 'rpointer.ice'
  character(*),parameter     :: u_FILE_u = &
       __FILE__

  type i2x_type
     real(r8), pointer ::  Si_imask(:)      => null()    
     real(r8), pointer ::  Si_ifrac(:)      => null()    
     real(r8), pointer ::  Si_t(:)          => null()        
     real(r8), pointer ::  Si_tref(:)       => null()     
     real(r8), pointer ::  Si_qref(:)       => null()     
     real(r8), pointer ::  Si_avsdr(:)      => null()    
     real(r8), pointer ::  Si_anidr(:)      => null()    
     real(r8), pointer ::  Si_avsdf(:)      => null()    
     real(r8), pointer ::  Si_anidf(:)      => null()    
     real(r8), pointer ::  Faii_swnet(:)    => null()  
     real(r8), pointer ::  Faii_sen(:)      => null()    
     real(r8), pointer ::  Faii_lat(:)      => null()    
     real(r8), pointer ::  Faii_lwup(:)     => null()   
     real(r8), pointer ::  Faii_evap(:)     => null()   
     real(r8), pointer ::  Faii_taux(:)     => null()   
     real(r8), pointer ::  Faii_tauy(:)     => null()   
     real(r8), pointer ::  Fioi_melth(:)    => null()  
     real(r8), pointer ::  Fioi_meltw(:)    => null()  
     real(r8), pointer ::  Fioi_swpen(:)    => null()  
     real(r8), pointer ::  Fioi_taux(:)     => null()   
     real(r8), pointer ::  Fioi_tauy(:)     => null()   
     real(r8), pointer ::  Fioi_salt(:)     => null()   
     real(r8), pointer ::  Fioi_bcpho(:)    => null()  
     real(r8), pointer ::  Fioi_bcphi(:)    => null()  
     real(r8), pointer ::  Fioi_flxdst(:)   => null() 
     real(r8), pointer ::  Si_ifrac_01(:,:) => null()
     real(r8), pointer ::  Fioi_swpen_ifrac_01(:,:) => null()  
  end type i2x_type
  type(i2x_type) :: i2x

  type x2i_type
     real(r8), pointer :: Faxa_swvdr(:)    => null()
     real(r8), pointer :: Faxa_swvdf(:)    => null()
     real(r8), pointer :: Faxa_swndr(:)    => null()
     real(r8), pointer :: Faxa_swndf(:)    => null()
     real(r8), pointer :: Fioo_q(:)        => null()
     real(r8), pointer :: Sa_z(:)          => null()
     real(r8), pointer :: Sa_u(:)          => null()
     real(r8), pointer :: Sa_v(:)          => null()
     real(r8), pointer :: Sa_ptem(:)       => null()
     real(r8), pointer :: Sa_shum(:)       => null()
     real(r8), pointer :: Sa_dens(:)       => null()
     real(r8), pointer :: Sa_tbot(:)       => null()
     real(r8), pointer :: So_s(:)          => null()
     real(r8), pointer :: Faxa_bcph(:,:)   => null()
     real(r8), pointer :: Faxa_ocph(:,:)   => null()
     real(r8), pointer :: Faxa_dstdry(:,:) => null()
     real(r8), pointer :: Faxa_dstwet(:,:) => null()
  end type x2i_type
  type(x2i_type) :: x2i

  logical :: flds_i2o_per_cat_mod

!===============================================================================
contains
!===============================================================================

  subroutine dice_comp_advertise(importState, exportState, flds_scalar_name, &
       flds_i2o_per_cat, fldsFrIce_num, fldsFrIce, fldsToIce_num, fldsToIce, rc)

    ! --------------------------------------------------------------
    ! advertitse import and export fields to mediator
    ! --------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    character(len=*)     , intent(in)    :: flds_scalar_name
    logical              , intent(in)    :: flds_i2o_per_cat
    integer              , intent(out)   :: fldsToIce_num
    integer              , intent(out)   :: fldsFrIce_num
    type (fld_list_type) , intent(out)   :: fldsToIce(:)
    type (fld_list_type) , intent(out)   :: fldsFrIce(:)
    integer              , intent(out)   :: rc

    ! local variables
    integer         :: n
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (trim(datamode) == 'NULL') then
       RETURN
    end if

    ! Advertise export fields

    fldsFrIce_num=1
    fldsFrIce(1)%stdname = trim(flds_scalar_name)

    call dshr_fld_add('Si_ifrac'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_imask'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_t'        , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_tref'     , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_qref'     , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_avsdr'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_anidr'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_avsdf'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Si_anidf'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_swnet'  , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_sen'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_lat'    , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_lwup'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_evap'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_taux'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Faii_tauy'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_melth'  , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_meltw'  , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_swpen'  , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_taux'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_tauy'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_salt'   , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_bcpho'  , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_bcphi'  , fldsFrIce_num, fldsFrIce)
    call dshr_fld_add('Fioi_flxdst' , fldsFrIce_num, fldsFrIce)
    if (flds_i2o_per_cat) then
       call dshr_fld_add('Si_ifrac_n'        , fldsFrIce_num, fldsFrIce, ungridded_lbound=1, ungridded_ubound=1)
       call dshr_fld_add('Fioi_swpen_ifrac_n', fldsFrIce_num, fldsFrIce, ungridded_lbound=1, ungridded_ubound=1)
    end if

    do n = 1,fldsFrIce_num
       call NUOPC_Advertise(exportState, standardName=fldsFrIce(n)%stdname, TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    ! Advertise import fields

    if (trim(datamode) /= 'NULL') then
       fldsToIce_num=1
       fldsToIce(1)%stdname = trim(flds_scalar_name)

       call dshr_fld_add('Faxa_swvdr' , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Faxa_swvdf' , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Faxa_swndr' , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Faxa_swndf' , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Fioo_q'     , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_z'       , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_u'       , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_v'       , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_ptem'    , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_shum'    , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_dens'    , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Sa_tbot'    , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('So_s'       , fldlist_num=fldsToIce_num, fldlist=fldsToIce)
       call dshr_fld_add('Faxa_bcph'  , fldlist_num=fldsToIce_num, fldlist=fldsToIce, ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fld_add('Faxa_ocph'  , fldlist_num=fldsToIce_num, fldlist=fldsToIce, ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fld_add('Faxa_dstdry', fldlist_num=fldsToIce_num, fldlist=fldsToIce, ungridded_lbound=1, ungridded_ubound=4)
       call dshr_fld_add('Faxa_dstwet', fldlist_num=fldsToIce_num, fldlist=fldsToIce, ungridded_lbound=1, ungridded_ubound=4)

       do n = 1,fldsToIce_num
          call NUOPC_Advertise(importState, standardName=fldsToIce(n)%stdname, TransferOfferGeomObject='will provide', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       enddo
    end if

    flds_i2o_per_cat_mod = flds_i2o_per_cat

  end subroutine dice_comp_advertise

  !===============================================================================

  subroutine dice_comp_setptrs(importState, exportState, rc)

    !-------------------
    ! Initialize contents of i2x and x2i
    !-------------------

    ! input/output parameters
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_StateItem_Flag) :: itemFlag
    integer                   :: n, lsize, kf
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Set pointers to exportState fields
    call state_getfldptr(exportState, fldname='Si_ifrac'    , fldptr1=i2x%Si_ifrac    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_imask'    , fldptr1=i2x%Si_imask    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_t'        , fldptr1=i2x%Si_t        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_tref'     , fldptr1=i2x%Si_tref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_qref'     , fldptr1=i2x%Si_qref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_avsdr'    , fldptr1=i2x%Si_avsdr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_anidr'    , fldptr1=i2x%Si_anidr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_avsdf'    , fldptr1=i2x%Si_avsdf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_anidf'    , fldptr1=i2x%Si_anidf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_swnet'  , fldptr1=i2x%Faii_swnet  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_sen'    , fldptr1=i2x%Faii_sen    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_lat'    , fldptr1=i2x%Faii_lat    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_lwup'   , fldptr1=i2x%Faii_lwup   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_evap'   , fldptr1=i2x%Faii_evap   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_taux'   , fldptr1=i2x%Faii_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_tauy'   , fldptr1=i2x%Faii_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_melth'  , fldptr1=i2x%Fioi_melth  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_meltw'  , fldptr1=i2x%Fioi_meltw  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_swpen'  , fldptr1=i2x%Fioi_swpen  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_taux'   , fldptr1=i2x%Fioi_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_tauy'   , fldptr1=i2x%Fioi_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_salt'   , fldptr1=i2x%Fioi_salt   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_bcpho'  , fldptr1=i2x%Fioi_bcpho  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_bcphi'  , fldptr1=i2x%Fioi_bcphi  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_flxdst' , fldptr1=i2x%Fioi_flxdst , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_i2o_per_cat_mod) then
       call state_getfldptr(exportState, fldname='Si_ifrac_01'        , fldptr2=i2x%Si_ifrac_01, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(exportState, fldname='Fioi_swpen_ifrac_01', fldptr2=i2x%Fioi_swpen_ifrac_01, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Set pointers to importState fields

    if (trim(datamode) /= 'NULL') then
       call state_getfldptr(importState, fldname='Faxa_swvdr'  , fldptr1=x2i%Faxa_swvdr  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_swvdf'  , fldptr1=x2i%Faxa_swvdf  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_swndr'  , fldptr1=x2i%Faxa_swndr  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_swndf'  , fldptr1=x2i%Faxa_swndf  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_bcph'   , fldptr2=x2i%Faxa_bcph   , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_ocph'   , fldptr2=x2i%Faxa_ocph   , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_dstdry' , fldptr2=x2i%Faxa_dstdry , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Faxa_dstwet' , fldptr2=x2i%Faxa_dstwet , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Fioo_q'      , fldptr1=x2i%Fioo_q      , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_z'        , fldptr1=x2i%Sa_z        , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_u'        , fldptr1=x2i%Sa_u        , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_v'        , fldptr1=x2i%Sa_v        , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_ptem'     , fldptr1=x2i%Sa_ptem     , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_tbot'     , fldptr1=x2i%Sa_tbot     , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_shum'     , fldptr1=x2i%Sa_shum     , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='Sa_dens'     , fldptr1=x2i%Sa_dens     , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, fldname='So_s'        , fldptr1=x2i%So_s        , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! initialize x2i arrays
    ! On initial call, x2i is unset, so set for the first use in generating the export state
    ! These values should have no impact on the solution!!

    x2i%Faxa_swvdr(:)    = 0._r8
    x2i%Faxa_swvdf(:)    = 0._r8
    x2i%Faxa_swndr(:)    = 0._r8
    x2i%Faxa_swndf(:)    = 0._r8
    x2i%Faxa_bcph(:,:)   = 0._r8
    x2i%Faxa_ocph(:,:)   = 0._r8
    x2i%Faxa_dstdry(:,:) = 0._r8
    x2i%Faxa_dstwet(:,:) = 0._r8
    x2i%Fioo_q(:)        = 0._r8
    x2i%Sa_z(:)          = 10.0_r8
    x2i%Sa_u(:)          = 5.0_r8
    x2i%Sa_v(:)          = 5.0_r8
    x2i%Sa_ptem(:)       = 260.0_r8
    x2i%Sa_tbot(:)       = 260.0_r8
    x2i%Sa_shum(:)       = 0.0014_r8
    x2i%Sa_dens(:)       = 1.3_r8
    x2i%So_s(:)          = 0._r8

  end subroutine dice_comp_setptrs

  !===============================================================================

  subroutine dice_comp_copy_streams()

    !-------------------
    ! Loop over streams, and map the correct stream field to the appropriate field pointer
    !-------------------

    ! local variables
    integer :: n, kf , nfld
    !-------------------------------------------------------------------------------

    do nfld = 1, SDICE%nstreams
       kf = mct_aVect_indexRA(SDICE%avs(nfld), 'ifrac', perrWith='quiet')
       if (kf /= 0) then
          do n = 1,size(i2x%Si_ifrac)
             i2x%Si_ifrac(n) = SDICE%avs(nfld)%rattr(kf,n)
          end do
       end if
    end do

  end subroutine dice_comp_copy_streams

  !===============================================================================

  subroutine dice_comp_init(mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon, calendar, mesh, nxg, nyg, rc)

    ! ----------------------------
    ! initialize dice model
    ! ----------------------------

    ! input/output parameters:
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: compid           ! mct comp id
    integer                , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name        ! fullname of current instance (ie. "lnd_0001")
    integer                , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: scmMode          ! single column mode
    real(R8)               , intent(in)    :: scmLat           ! single column lat
    real(R8)               , intent(in)    :: scmLon           ! single column lon
    character(len=*)       , intent(in)    :: calendar         ! calendar type
    type(ESMF_Mesh)        , intent(in)    :: mesh             ! ESMF dice mesh
    integer                , intent(out)   :: nxg, nyg
    integer                , intent(out)   :: rc

    !--- local variables ---
    integer                      :: n,k            ! generic counters
    integer                      :: ierr           ! error code
    integer                      :: lsize          ! local size
    integer                      :: kfld           ! field reference
    logical                      :: exists,exists1 ! file existance logical
    integer                      :: nu             ! unit number
    type(ESMF_DistGrid)          :: distGrid
    integer, allocatable, target :: gindex(:)
    integer                      :: dimCount
    integer                      :: tileCount
    integer                      :: deCount
    integer                      :: gsize
    integer, allocatable         :: elementCountPTile(:)
    integer, allocatable         :: indexCountPDE(:,:)
    integer                      :: spatialDim
    integer                      :: numOwnedElements
    real(R8), pointer            :: ownedElemCoords(:)
    character(*), parameter      :: F00   = "('(dice_comp_init) ',8a)"
    character(*), parameter      :: F01   = "('(dice_comp_init) ',a,2f10.4)"
    character(*), parameter      :: subName = "(dice_comp_init) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------------------------------------------------------------
    ! Initialize PIO
    !----------------------------------------------------------------------------

    call t_startf('DICE_INIT')

    call shr_strdata_pioinit(SDICE, compid)

    !----------------------------------------------------------------------------
    ! Create a data model global segmap
    !----------------------------------------------------------------------------

    call t_startf('dice_strdata_init')

    if (my_task == master_task) write(logunit,F00) ' initialize SDICE gsmap'

    ! obtain the distgrid from the mesh that was read in
    call ESMF_MeshGet(Mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determin local size on my processor
    call ESMF_distGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global index space for my processor
    allocate(gindex(lsize))
    call ESMF_distGridGet(distGrid, localDe=0, seqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global size of distgrid
    call ESMF_distGridGet(distGrid, dimCount=dimCount, deCount=deCount, tileCount=tileCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    gsize = 0
    do n = 1,size(elementCountPTile)
       gsize = gsize + elementCountPTile(n)
    end do
    deallocate(elementCountPTile)

    ! create the data model gsmap given the local size, global size and gindex
    call mct_gsMap_init( SDICE%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    !----------------------------------------------------------------------------
    ! Initialize SDICE
    !----------------------------------------------------------------------------

    ! The call to shr_strdata_init_model_domain creates the SDICE%gsmap which
    ! is a '2d1d' decommp (1d decomp of 2d grid) and also create SDICE%grid

    SDICE%calendar = trim(shr_cal_calendarName(trim(calendar)))

    if (scmmode) then
       if (my_task == master_task) write(logunit,F01) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init_model_domain(SDICE, mpicom, compid, my_task, &
            scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, gsmap=SDICE%gsmap)
    else
       call shr_strdata_init_model_domain(SDICE, mpicom, compid, my_task, gsmap=SDICE%gsmap)
    end if

    if (my_task == master_task) then
       call shr_strdata_print(SDICE,'SDICE data')
    endif

    ! obtain mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(xc(numOwnedElements), yc(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (numOwnedElements /= lsize) then
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if
    do n = 1,lsize
       xc(n) = ownedElemCoords(2*n-1)
       yc(n) = ownedElemCoords(2*n)
    end do

    ! error check that mesh lats and lons correspond to those on the input domain file
    index_lon = mct_aVect_indexRA(SDICE%grid%data,'lon')
    do n = 1, lsize
       if (abs(mod(SDICE%grid%data%rattr(index_lon,n) - xc(n),360.0_R8)) > 1.e-4) then
          write(6,*)'ERROR: lon diff = ',abs(SDICE%grid%data%rattr(index_lon,n) -  xc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDICE%grid%data%rattr(index_lon,n) = xc(n) ! overwrite ggrid with mesh data
    end do
    index_lat = mct_aVect_indexRA(SDICE%grid%data,'lat')
    do n = 1, lsize
       if (abs( SDICE%grid%data%rattr(index_lat,n) -  yc(n)) > 1.e-4) then
          write(6,*)'ERROR: lat diff = ',abs(SDICE%grid%data%rattr(index_lat,n) -  yc(n)),' too large'
          call shr_sys_abort()
       end if
       !SDICE%grid%data%rattr(index_lat,n) = yc(n) ! overwrite ggrid with mesh data
    end do

    ! Note that the module array, imask, does not change after initialization
    kfld = mct_aVect_indexRA(SDICE%grid%data,'mask')
    i2x%Si_imask(:) = SDICE%grid%data%rAttr(kfld,:)

    allocate(imask(lsize))
    imask(:) = nint(SDICE%grid%data%rAttr(kfld,:))

    if (my_task == master_task) then
       call shr_strdata_print(SDICE,'SDICE data')
    endif

    call t_stopf('dice_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize SDICE attributes for streams and mapping of streams to model domain
    !----------------------------------------------------------------------------

    call shr_strdata_init_streams(SDICE, compid, mpicom, my_task)
    call shr_strdata_init_mapping(SDICE, compid, mpicom, my_task)

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dice_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'

    allocate(water(lsize))
    allocate(tfreeze(lsize))
    ! allocate(iFrac0(lsize))

    call t_stopf('dice_initmctavs')

    nxg = SDICE%nxg
    nyg = SDICE%nyg

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       exists = .false.
       exists1 = .false.
       if (trim(rest_file)      == trim(nullstr) .and. &
           trim(rest_file_strm) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer = ',trim(rpfile)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (exists) then
                open(newunit=nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
                read(nu,'(a)') rest_file
                read(nu,'(a)') rest_file_strm
                close(nu)
                inquire(file=trim(rest_file_strm),exist=exists)
                inquire(file=trim(rest_file),exist=exists1)
             endif
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm,mpicom,'rest_file_strm')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
       endif

       call shr_mpi_bcast(exists,mpicom,'exists')
       call shr_mpi_bcast(exists1,mpicom,'exists1')

       if (exists1) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
          call shr_pcdf_readwrite('read',SDICE%pio_subsystem, SDICE%io_type, &
               trim(rest_file), mpicom, gsmap=SDICE%gsmap, rf1=water, rf1n='water', io_format=SDICE%io_format)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file)
       endif

       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDICE,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
    endif

    call t_stopf('DICE_INIT')

  end subroutine dice_comp_init

  !===============================================================================

  subroutine dice_comp_run(mpicom, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       calendar, modeldt, target_ymd, target_tod, cosArg, case_name, rc )

    ! --------------------------
    ! run method for dice model
    ! --------------------------

    ! input/output parameters:
    integer                , intent(in)    :: mpicom               ! mpi communicator
    integer                , intent(in)    :: my_task              ! my task in mpi communicator mpicom
    integer                , intent(in)    :: master_task          ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix          ! char string associated with instance
    integer                , intent(in)    :: logunit              ! logging unit number
    logical                , intent(in)    :: read_restart         ! start from restart
    logical                , intent(in)    :: write_restart        ! restart now
    character(len=*)       , intent(in)    :: calendar
    integer                , intent(in)    :: modeldt
    integer                , intent(in)    :: target_ymd
    integer                , intent(in)    :: target_tod
    real(R8)               , intent(in)    :: cosarg               ! for setting ice temp pattern
    character(len=*)       , intent(in), optional :: case_name     ! case name
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n,nfld            ! indices
    integer           :: lsize             ! size of attr vect
    real(R8)          :: dt                ! timestep
    integer           :: nu                ! unit number
    real(R8)          :: qmeltall          ! q that would melt all accumulated water
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dice_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dice_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: F0D   = "('(dice_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--------------------
    ! ADVANCE ICE
    !--------------------

    call t_startf('DICE_RUN')
    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice')

    dt = modeldt * 1.0_r8

    !--------------------
    ! Advance streams
    !--------------------

    if (trim(datamode) /= 'NULL') then
       call t_startf('dice_strdata_advance')
       call shr_strdata_advance(SDICE, target_ymd, target_tod, mpicom, 'dice')
       call t_stopf('dice_strdata_advance')
    end if

    !--------------------
    ! Copy all fields from streams to i2x
    !--------------------

    if (trim(datamode) /= 'NULL') then
       call t_barrierf('dice_comp_copy_streams_BARRIER', mpicom)
       call t_startf('dice_comp_copy_streams')
       call dice_comp_copy_streams()
       call t_stopf('dice_comp_copy_streams')
    endif

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    lsize = size(i2x%Si_ifrac)

    call t_startf('dice_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       if (firstcall .and. .not. read_restart) then
          ! iFrac0 = iFrac  ! previous step's ice fraction
          water(:) = 0.0_R8 ! previous step's water accumulation
          where (i2x%Si_ifrac(:) > 0.0_R8) water(:) = flux_Qacc0
       endif

       tfreeze(:) = shr_frz_freezetemp(x2i%So_s(:)) + tFrz ! convert to Kelvin

       do n = 1,lsize

          !--- fix erroneous iFrac ---
          i2x%Si_ifrac(n) = min(1.0_R8,max(0.0_R8,i2x%Si_ifrac(n)))

          !--- fabricate ice surface T, fix erroneous iFrac ---
          if ( yc(n) > 0.0_R8) then
             i2x%Si_t(n) = 260.0_R8 + 10.0_R8*cos(cosArg)
          else
             i2x%Si_t(n) = 260.0_R8 - 10.0_R8*cos(cosArg)
          end if

          !--- set albedos (constant) ---
          i2x%Si_avsdr(n) = ax_vsdr
          i2x%Si_anidr(n) = ax_nidr
          i2x%Si_avsdf(n) = ax_vsdf
          i2x%Si_anidf(n) = ax_nidf

          !--- swnet is sent to cpl as a diagnostic quantity only ---
          !--- newly recv'd swdn goes with previously sent albedo ---
          !--- but albedos are (currently) time invariant         ---
          i2x%Faii_swnet(n)   = (1.0_R8 - i2x%Si_avsdr(n))*x2i%Faxa_swvdr(n) &
                              + (1.0_R8 - i2x%Si_anidr(n))*x2i%Faxa_swndr(n) &
                              + (1.0_R8 - i2x%Si_avsdf(n))*x2i%Faxa_swvdf(n) &
                              + (1.0_R8 - i2x%Si_anidf(n))*x2i%Faxa_swndf(n)

          !--- compute melt/freeze water balance, adjust iFrac  -------------
          if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
             i2x%Fioi_melth(n) = min(x2i%Fioo_q(n),0.0_R8 )          ! q<0 => melt potential
             i2x%Fioi_melth(n) = max(i2x%Fioi_melth(n),Flux_Qmin   ) ! limit the melt rate
             i2x%Fioi_meltw(n) =    -i2x%Fioi_melth(n)/latice        ! corresponding water flux

          else                                 ! Q accumulation option is ON
             !--------------------------------------------------------------
             ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
             ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
             ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
             ! 2b) Q>0 & iFrac = 0  =>  accumulated water
             !--------------------------------------------------------------

             if ( x2i%Fioo_q(n) <  0.0_R8 ) then ! Q<0 => melt
                if (i2x%Si_ifrac(n) > 0.0_R8 ) then
                   i2x%Fioi_melth(n) = i2x%Si_ifrac(n)*max(x2i%Fioo_q(n),Flux_Qmin)
                   i2x%Fioi_meltw(n) =    -i2x%Fioi_melth(n)/latice
                   !  water(n) = < don't change this value >
                else
                   Qmeltall   = -water(n)*latice/dt
                   i2x%Fioi_melth(n) = max(x2i%Fioo_q(n), Qmeltall, Flux_Qmin )
                   i2x%Fioi_meltw(n) = -i2x%Fioi_melth(n)/latice
                   water(n) =  water(n) - i2x%Fioi_meltw(n)*dt
                end if
             else                       ! Q>0 => freeze
                if (i2x%Si_ifrac(n) > 0.0_R8 ) then
                   i2x%Fioi_melth(n) = 0.0_R8
                   i2x%Fioi_meltw(n) = 0.0_R8
                   water(n) = 0.0_R8
                else
                   i2x%Fioi_melth(n) = 0.0_R8
                   i2x%Fioi_meltw(n) = 0.0_R8
                   water(n) = water(n) + dt*x2i%Fioo_q(n)/latice
                end if
             end if

             if (water(n) < 1.0e-16_R8 ) water(n) = 0.0_R8

             !--- non-zero water => non-zero iFrac ---
             if (i2x%Si_ifrac(n) <= 0.0_R8  .and.  water(n) > 0.0_R8) then
                i2x%Si_ifrac(n) = min(1.0_R8,water(n)/waterMax)
                ! i2x%Si_t(n) = tfreeze(n)     ! T can be above freezing?!?
             end if

             !--- cpl multiplies Fioi_melth & Fioi_meltw by iFrac ---
             !--- divide by iFrac here => fixed quantity flux (not per area) ---
             if (i2x%Si_ifrac(n) > 0.0_R8) then
                i2x%Si_ifrac(n) = max( 0.01_R8, i2x%Si_ifrac(n)) ! min iFrac
                i2x%Fioi_melth(n) = i2x%Fioi_melth(n)/i2x%Si_ifrac(n)
                i2x%Fioi_meltw(n) = i2x%Fioi_meltw(n)/i2x%Si_ifrac(n)
             else
                i2x%Fioi_melth(n) = 0.0_R8
                i2x%Fioi_meltw(n) = 0.0_R8
             end if
          end if

          !--- modify T wrt iFrac: (iFrac -> 0) => (T -> tfreeze) ---
          i2x%Si_t(n) = tfreeze(n) + i2x%Si_ifrac(n)*(i2x%Si_t(n)-tfreeze(n))
       end do

       ! compute ice/ice surface fluxes
       call dice_flux_atmice( &
            imask         ,x2i%Sa_z      ,x2i%Sa_u      ,x2i%Sa_v     , &
            x2i%Sa_ptem   ,x2i%Sa_shum   ,x2i%Sa_dens   ,x2i%Sa_tbot  , &
            i2x%Si_t      ,i2x%Faii_sen  ,i2x%Faii_lat  ,i2x%Faii_lwup, &
            i2x%Faii_evap ,i2x%Faii_taux ,i2x%Faii_tauy ,i2x%Si_tref  , &
            i2x%Si_qref   ,logunit )

       ! compute ice/oce surface fluxes (except melth & meltw, see above)
       do n=1,lsize
          if (imask(n) == 0) then
             i2x%Fioi_swpen(n) = spval
             i2x%Fioi_melth(n) = spval
             i2x%Fioi_meltw(n) = spval
             i2x%Fioi_salt (n) = spval
             i2x%Fioi_taux(n)  = spval
             i2x%Fioi_tauy(n)  = spval
             i2x%Si_ifrac(n)   = 0.0_R8
          else
             !--- penetrating short wave ---
             i2x%Fioi_swpen(n) = max(0.0_R8, flux_swpf*i2x%Faii_swnet(n) ) ! must be non-negative

             !--- i/o surface stress ( = atm/ice stress) ---
             i2x%Fioi_taux(n) = i2x%Faii_taux(n)
             i2x%Fioi_tauy(n) = i2x%Faii_tauy(n)

             !--- salt flux ---
             i2x%Fioi_salt(n) = 0.0_R8
          end if
          ! !--- save ifrac for next timestep
          ! iFrac0(n) = i2x%Si_ifrac(n)
       end do

       ! Compute outgoing aerosol fluxes
       do n = 1,lsize
          i2x%Fioi_bcpho(n) = x2i%Faxa_bcph(2,n)
          i2x%Fioi_bcphi(n) = x2i%Faxa_bcph(1,n) + x2i%Faxa_bcph(3,n)
          i2x%Fioi_flxdst(n) =  x2i%Faxa_dstdry(1,n) + x2i%Faxa_dstwet(1,n) &
                              + x2i%Faxa_dstdry(2,n) + x2i%Faxa_dstwet(2,n) &
                              + x2i%Faxa_dstdry(3,n) + x2i%Faxa_dstwet(3,n) &
                              + x2i%Faxa_dstdry(4,n) + x2i%Faxa_dstwet(4,n)
       end do

    end select

    !-------------------------------------------------
    ! optional per thickness category fields
    !-------------------------------------------------

    if (flds_i2o_per_cat_mod) then
       do n=1,lsize
          i2x%Si_iFrac_01(1,n)         = i2x%Si_ifrac(n)
          i2x%Fioi_swpen_iFrac_01(1,n) = i2x%Fioi_swpen(n) * i2x%Si_ifrac(n)
       end do
    end if

    call t_stopf('dice_datamode')

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('dice_restart')
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
       write(rest_file,"(6a)") &
            trim(case_name), '.dice',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.dice',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          open(newunit=nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
       endif
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),target_ymd,target_tod
       call shr_pcdf_readwrite('write',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file), mpicom, SDICE%gsmap, clobber=.true., rf1=water, rf1n='water')
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm),SDICE,mpicom,trim(case_name),'SDICE strdata')
       call t_stopf('dice_restart')
    endif

    call t_stopf('dice')

    firstcall = .false.

    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

end module dice_comp_mod
