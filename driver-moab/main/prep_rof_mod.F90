module prep_rof_mod

#include "shr_assert.h"
  use shr_kind_mod,     only: R8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: cxx => SHR_KIND_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_lnd, num_inst_rof, num_inst_frc, num_inst_atm, num_inst_ocn
  use seq_comm_mct,     only: CPLID, ROFID, logunit
  use seq_comm_mct,     only: mblxid   ! iMOAB id for land on coupler (read now from h5m file)
  use seq_comm_mct,     only: mbaxid   ! iMOAB id for atm migrated mesh to coupler pes (migrate either mhid or mhpgx, depending on atm_pg_active)
  use seq_comm_mct,     only: mbrxid   ! iMOAB id of moab rof read on couple pes
  use seq_comm_mct,     only: mbintxar ! iMOAB id for intx mesh between atm and river
  use seq_comm_mct,     only: mboxid   
  use seq_comm_mct,     only: mbintxlr ! iMOAB id for intx mesh between land and river
  use seq_comm_mct,     only : atm_pg_active  ! whether the atm uses FV mesh or not ; made true if fv_nphys > 0
  !use dimensions_mod,   only : np     ! for atmosphere degree 
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
  use component_type_mod, only: rof, lnd, atm, ocn
  use component_type_mod, only: ocn ! used for context for projection towards ocean from rof !
  use prep_lnd_mod, only: prep_lnd_get_mapper_Fr2l
  use map_lnd2rof_irrig_mod, only: map_lnd2rof_irrig
  use seq_comm_mct,     only: mb_rof_aream_computed  ! signal

  use iso_c_binding
#ifdef MOABCOMP
  use component_type_mod, only:  compare_mct_av_moab_tag
#endif

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_rof_init
  public :: prep_rof_mrg

#ifdef HAVE_MOAB
  public :: prep_rof_mrg_moab
  public :: prep_rof_accum_lnd_moab  
  public :: prep_rof_accum_atm_moab
  public :: prep_rof_accum_ocn_moab
  public :: prep_rof_accum_avg_moab
#endif

  public :: prep_rof_accum_lnd
  public :: prep_rof_accum_atm
  public :: prep_rof_accum_ocn
  public :: prep_rof_accum_avg

  public :: prep_rof_calc_l2r_rx
  public :: prep_rof_calc_a2r_rx
  public :: prep_rof_calc_o2r_rx

  public :: prep_rof_get_l2racc_lx
  public :: prep_rof_get_l2racc_lx_cnt
  public :: prep_rof_get_o2racc_ox
  public :: prep_rof_get_o2racc_ox_cnt
  public :: prep_rof_get_mapper_Fl2r
  public :: prep_rof_get_a2racc_ax
  public :: prep_rof_get_a2racc_ax_cnt
  public :: prep_rof_get_mapper_Sa2r
  public :: prep_rof_get_mapper_Fa2r

  public :: prep_rof_get_sharedFieldsOcnRof
  public :: prep_rof_get_o2racc_om ! return a pointer to private array !!! 
  public :: prep_rof_get_o2racc_om_cnt

  public :: prep_rof_get_l2racc_lm_cnt
  public :: prep_rof_get_l2racc_lm
  public :: prep_rof_get_sharedFieldsLndRof
  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_rof_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sa2r
  type(seq_map), pointer :: mapper_Fa2r
  type(seq_map), pointer :: mapper_Fl2r
  type(seq_map), pointer :: mapper_So2r

  ! attribute vectors
  type(mct_aVect), pointer :: l2r_rx(:)
  type(mct_aVect), pointer :: a2r_rx(:)
  type(mct_aVect), pointer :: o2r_rx(:)

  ! accumulation variables
  type(mct_aVect), pointer :: l2racc_lx(:)   ! lnd export, lnd grid, cpl pes
  integer        , target  :: l2racc_lx_cnt  ! l2racc_lx: number of time samples accumulated
  type(mct_aVect), pointer :: a2racc_ax(:)   ! atm export, atm grid, cpl pes
  integer        , target  :: a2racc_ax_cnt  ! a2racc_ax: number of time samples accumulated 
  type(mct_aVect), pointer :: o2racc_ox(:)   ! ocn export, ocn grid, cpl pes
  integer        , target  :: o2racc_ox_cnt  ! o2racc_ox: number of time samples accumulated

  ! accumulation variables over moab fields 
  character(CXX)                        :: sharedFieldsLndRof ! used in moab to define l2racc_lm
  real (kind=R8) , allocatable, private, target :: l2racc_lm(:,:)   ! lnd export, lnd grid, cpl pes
  real (kind=R8) , allocatable, private :: l2x_lm2(:,:)  ! basically l2x_lm, but in another copy, on rof module
  integer        , target  :: l2racc_lm_cnt  ! l2racc_lm: number of time samples accumulated
  integer :: nfields_sh_lr ! number of fields in sharedFieldsLndRof
  integer :: lsize_lm ! size of land in moab, local

  character(CXX)       :: sharedFieldsAtmRof ! used in moab to define a2racc_am
  real (kind=R8) , allocatable, private ::  a2racc_am(:,:)   ! atm export, atm grid, cpl pes
  real (kind=R8) , allocatable, private :: a2x_am2(:,:)  ! basically a2x_am, but in another copy, on rof module
  integer        , target  :: a2racc_am_cnt  ! a2racc_am: number of time samples accumulated 
  integer :: nfields_sh_ar ! number of fields in sharedFieldsAtmRof
  integer :: lsize_am ! size of atm in moab, local

  character(CXX)       :: sharedFieldsOcnRof ! used in moab to define o2racc_om
  real (kind=R8) , allocatable, private, target ::  o2racc_om(:,:)   ! ocn export, ocn grid, cpl pes
  real (kind=R8) , allocatable, private :: o2r_om2(:,:)  ! basically o2x_om, but in another copy, on rof module, only shared with rof
  integer        , target  :: o2racc_om_cnt  ! o2racc_om: number of time samples accumulated
  integer :: nfields_sh_or ! number of fields in sharedFieldsOcnRof
  integer :: lsize_om ! size of ocn in moab, local

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator

  ! field names and lists, for fields that need to be treated specially
  character(len=*), parameter :: irrig_flux_field = 'Flrl_irrig'
  ! fluxes mapped from lnd to rof that don't need any special handling
  character(CXX) :: lnd2rof_normal_fluxes
  ! whether the model is being run with a separate irrigation field
  logical :: have_irrig_field
  ! samegrid atm and lnd
  logical :: samegrid_al   ! samegrid atm and lnd

  ! moab stuff
   real (kind=R8) , allocatable, private :: fractions_rm (:,:) ! will retrieve the fractions from rof, and use them
  !  they were init with 
  ! character(*),parameter :: fraclist_r = 'lfrac:lfrin:rfrac'  in moab, on the fractions 
  real (kind=R8) , allocatable, private :: x2r_rm (:,:) ! result of merge
  real (kind=R8) , allocatable, private :: a2x_rm (:,:)
  real (kind=R8) , allocatable, private :: l2x_rm (:,:)

  !================================================================================================

contains

  !================================================================================================

  subroutine prep_rof_init(infodata, lnd_c2_rof, atm_c2_rof, ocn_c2_rof)

   use iMOAB, only: iMOAB_ComputeMeshIntersectionOnSphere, iMOAB_RegisterApplication, &
      iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_ComputeCommGraph, & 
      iMOAB_ComputeScalarProjectionWeights, iMOAB_GetMeshInfo
    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: lnd_c2_rof ! .true.  => lnd to rof coupling on
    logical                 , intent(in)    :: atm_c2_rof ! .true.  => atm to rof coupling on
    logical                 , intent(in)    :: ocn_c2_rof ! .true.  => ocn to rof coupling on
    !
    ! Local Variables
    integer                     :: lsize_r
    integer                     :: lsize_l
    integer                     :: lsize_a
    integer                     :: lsize_o
    integer                     :: eli, eri, eai, eoi
    logical                     :: samegrid_lr   ! samegrid lnd and rof
    logical                     :: samegrid_ar   ! samegrid atm and rof
    logical                     :: samegrid_ro   ! samegrid ocn and rof
    logical                     :: esmf_map_flag ! .true. => use esmf for mapping
    logical                     :: rof_present   ! .true.  => rof is present
    logical                     :: lnd_present   ! .true.  => lnd is present
    logical                     :: atm_present   ! .true.  => atm is present
    logical                     :: ocn_present   ! .true.  => ocn is present
    logical                     :: iamroot_CPLID ! .true. => CPLID masterproc
    character(CL)               :: atm_gnam      ! atm grid
    character(CL)               :: lnd_gnam      ! lnd grid
    character(CL)               :: rof_gnam      ! rof grid
    character(CL)               :: ocn_gnam      ! ocn grid
    type(mct_aVect) , pointer   :: l2x_lx
    type(mct_aVect) , pointer   :: a2x_ax
    type(mct_aVect) , pointer   :: o2x_ox
    type(mct_aVect) , pointer   :: x2r_rx
    integer                     :: index_irrig
    character(*)    , parameter :: subname = '(prep_rof_init)'
    character(*)    , parameter :: F00 = "('"//subname//" : ', 4A )"

    ! MOAB stuff
   integer                  :: ierr, idintx, rank
   character*32             :: appname, outfile, wopts, lnum
   character*32             :: dm1, dm2, dofnameS, dofnameT, wgtIdef
   integer                  :: orderS, orderT, volumetric, noConserve, validate, fInverseDistanceMap
   integer                  :: fNoBubble, monotonicity
! will do comm graph over coupler PES, in 2-hop strategy
   integer                  :: mpigrp_CPLID ! coupler pes group, used for comm graph phys <-> atm-ocn

   integer                  :: type1, type2 ! type for computing graph; should be the same type for ocean, 3 (FV)
   integer                  :: tagtype, numco, tagindex
   character(CXX)           :: tagName
   integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info

    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         rof_present=rof_present       , &
         lnd_present=lnd_present       , &
         atm_present=atm_present       , &
         ocn_present=ocn_present       , &
         lnd_gnam=lnd_gnam             , &
         atm_gnam=atm_gnam             , &
         rof_gnam=rof_gnam             )

    allocate(mapper_Sa2r)
    allocate(mapper_Fa2r)
    allocate(mapper_Fl2r)
    allocate(mapper_So2r)

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
#ifdef HAVE_MOAB
       ! this l2racc_lm will be over land size ? 
       sharedFieldsLndRof=''
       nfields_sh_lr = mct_aVect_nRAttr(l2racc_lx(1))
       if( nfields_sh_lr /= 0 ) sharedFieldsLndRof=trim( mct_aVect_exportRList2c(l2racc_lx(1)) )
       tagname = trim(sharedFieldsLndRof)//C_NULL_CHAR
       ! find the size of land mesh locally
       ! find out the number of local elements in moab mesh lnd instance on coupler
       ierr  = iMOAB_GetMeshInfo ( mblxid, nvert, nvise, nbl, nsurf, nvisBC )
       if (ierr .ne. 0) then
          write(logunit,*) subname,' error in getting info '
          call shr_sys_abort(subname//' error in getting info ')
       endif
       ! land is fully cell now
       lsize_lm = nvise(1)
       if(iamroot_CPLID) then
          write(logunit,*) subname,' number of fields=', nfields_sh_lr
          write(logunit,*) subname,' sharedFieldsLndRof=', trim(sharedFieldsLndRof)      
          write(logunit,*) subname,'  seq_flds_l2x_fluxes_to_rof=', trim(seq_flds_l2x_fluxes_to_rof)
          write(logunit,*) subname,' lsize_lm=', lsize_lm
       endif
       allocate(l2racc_lm(lsize_lm, nfields_sh_lr))
       allocate(l2x_lm2(lsize_lm, nfields_sh_lr)) ! this will be obtained from land instance
       l2racc_lm(:,:) = 0.
       l2racc_lm_cnt = 0
#endif
       allocate(l2r_rx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(l2r_rx(eri), rList=seq_flds_l2x_fluxes_to_rof, lsize=lsize_r)
          call mct_avect_zero(l2r_rx(eri))
       end do

       samegrid_lr = .true.
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.
       samegrid_al = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.

       if (lnd_c2_rof) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fl2r'
          end if
          call seq_map_init_rcfile(mapper_Fl2r, lnd(1), rof(1), &
               'seq_maps.rc','lnd2rof_fmapname:','lnd2rof_fmaptype:',samegrid_lr, &
               string='mapper_Fl2r initialization', esmf_map=esmf_map_flag)
! similar to a2r, from below 
#ifdef HAVE_MOAB
          ! Call moab intx only if land and river are init in moab
          if ((mblxid .ge. 0) .and.  (mbrxid .ge. 0)) then
            appname = "LND_ROF_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between lnd and rof mesh
            idintx = 100*lnd(1)%cplcompid + rof(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxlr)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering lnd rof intx'
              call shr_sys_abort(subname//' ERROR in registering lnd rof intx')
            endif
            tagname = trim(seq_flds_l2x_fluxes_to_rof)//C_NULL_CHAR
            tagtype = 1 ! dense
            numco = 1 ! 
            ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags for seq_flds_a2x_fields on rof cpl'
               call shr_sys_abort(subname//' ERROR in  defining tags for seq_flds_a2x_fields on rof cpl')
            endif
            call seq_comm_getData(CPLID ,mpigrp=mpigrp_CPLID)
            if (samegrid_lr) then
               ! the same mesh , lnd and rof use the same dofs, but restricted 
               ! we do not compute intersection, so we will have to just send data from lnd to rof and viceversa, by GLOBAL_ID matching
               ! so we compute just a comm graph, between lnd and rof dofs, on the coupler; target is rof 
               ! land is full mesh
               type1 = 3; !  full mesh for land now
               type2 = 3;  ! fv for target rof
               ierr = iMOAB_ComputeCommGraph( mblxid, mbrxid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        lnd(1)%cplcompid, rof(1)%cplcompid)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing comm graph , lnd-rof'
                  call shr_sys_abort(subname//' ERROR in computing comm graph , lnd-rof')
               endif
               ! context for rearrange is target in this case
               if ( mapper_Fl2r%src_mbid .gt. -1 ) then
                  if (iamroot_CPLID) then
                        write(logunit,F00) 'overwriting '//trim(mapper_Fl2r%mbname) &
                              //' mapper_Fl2r'
                  endif
               endif
               mapper_Fl2r%src_mbid = mblxid
               mapper_Fl2r%tgt_mbid = mbrxid
               mapper_Fl2r%src_context = lnd(1)%cplcompid
               mapper_Fl2r%intx_context = rof(1)%cplcompid
            else
               ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mblxid, mbrxid, mbintxlr)
               if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing  land rof intx'
               call shr_sys_abort(subname//' ERROR in computing land rof intx')
               endif
               if (iamroot_CPLID) then
               write(logunit,*) 'iMOAB intersection between  land and rof with id:', idintx
               end if
               ! we also need to compute the comm graph for the second hop, from the lnd on coupler to the 
               ! lnd for the intx lnd-rof context (coverage)
               !    
               type1 = 3 ! land is FV now on coupler side
               type2 = 3;

               ierr = iMOAB_ComputeCommGraph( mblxid, mbintxlr, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                          lnd(1)%cplcompid, idintx)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing comm graph for second hop, lnd-rof'
                  call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, lnd-rof')
               endif
               ! now take care of the mapper 
               if ( mapper_Fl2r%src_mbid .gt. -1 ) then
                  if (iamroot_CPLID) then
                        write(logunit,F00) 'overwriting '//trim(mapper_Fl2r%mbname) &
                              //' mapper_Fl2r'
                  endif
               endif
               mapper_Fl2r%src_mbid = mblxid
               mapper_Fl2r%tgt_mbid = mbrxid
               mapper_Fl2r%intx_mbid = mbintxlr 
               mapper_Fl2r%src_context = lnd(1)%cplcompid
               mapper_Fl2r%intx_context = idintx
               wgtIdef = 'scalar'//C_NULL_CHAR
               mapper_Fl2r%weight_identifier = wgtIdef
               mapper_Fl2r%mbname = 'mapper_Fl2r'
               ! because we will project fields from lnd to rof grid, we need to define 
               !  the l2x fields to rof grid on coupler side
               
               volumetric = 0 ! can be 1 only for FV->DGLL or FV->CGLL; 
               
               dm1 = "fv"//C_NULL_CHAR
               dofnameS="GLOBAL_ID"//C_NULL_CHAR
               orderS = 1 !  fv-fv
            
               dm2 = "fv"//C_NULL_CHAR
               dofnameT="GLOBAL_ID"//C_NULL_CHAR
               orderT = 1  !  not much arguing
               fNoBubble = 1
               monotonicity = 0 !
               noConserve = 0
               validate = 0
               fInverseDistanceMap = 0
               if (iamroot_CPLID) then
                  write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxlr=', mbintxlr, ' wgtIdef=', wgtIdef, &
                     'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT)
               endif
               ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxlr, wgtIdef, &
                                                trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT) )
               mb_rof_aream_computed = .true. ! signal
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing lr weights '
                  call shr_sys_abort(subname//' ERROR in computing lr weights ')
               endif

#ifdef MOABDEBUG
               wopts = C_NULL_CHAR
               call shr_mpi_commrank( mpicom_CPLID, rank )
               if (rank .lt. 5) then
               write(lnum,"(I0.2)")rank !
               outfile = 'intx_lr_'//trim(lnum)// '.h5m' // C_NULL_CHAR
               ierr = iMOAB_WriteMesh(mbintxlr, outfile, wopts) ! write local intx file
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in writing intx lr file '
                  call shr_sys_abort(subname//' ERROR in writing intx lr file ')
               endif
               endif
#endif
            end if ! if ((mblxid .ge. 0) .and.  (mbrxid .ge. 0))
         endif ! samegrid_lr
#endif
          ! We'll map irrigation specially, so exclude this from the list of l2r fields
          ! that are mapped "normally".
          !
          ! (This listDiff works even if have_irrig_field is false.)
          call shr_string_listDiff( &
               list1 = seq_flds_l2x_fluxes_to_rof, &
               list2 = irrig_flux_field, &
               listout = lnd2rof_normal_fluxes)
       endif ! if (lnd_c2_rof) then
       call shr_sys_flush(logunit)

    end if ! if (rof_present .and. lnd_present) then

    if (rof_present .and. atm_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       lsize_r = mct_aVect_lsize(x2r_rx)

       a2x_ax => component_get_c2x_cx(atm(1))
       lsize_a = mct_aVect_lsize(a2x_ax)

       allocate(a2racc_ax(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_initSharedFields(a2x_ax, x2r_rx, a2racc_ax(eai), lsize=lsize_a)
          call mct_aVect_zero(a2racc_ax(eai))
       end do
       a2racc_ax_cnt = 0
#ifdef HAVE_MOAB
       ! this a2racc_am will be over atm size 
       sharedFieldsAtmRof=''
       nfields_sh_ar = mct_aVect_nRAttr(a2racc_ax(1))
       if (nfields_sh_ar /= 0 ) sharedFieldsAtmRof = trim( mct_aVect_exportRList2c(a2racc_ax(1)) )
       tagname = trim(sharedFieldsAtmRof)//C_NULL_CHAR
       ! find the size of atm mesh locally
       ! find out the number of local elements in moab mesh atm instance on coupler
       ierr  = iMOAB_GetMeshInfo ( mbaxid, nvert, nvise, nbl, nsurf, nvisBC )
       if (ierr .ne. 0) then
          write(logunit,*) subname,' error in getting info '
          call shr_sys_abort(subname//' error in getting info ')
       endif
       ! land is fully cell now
       lsize_am = nvise(1)
       allocate(a2racc_am(lsize_am, nfields_sh_ar))
       allocate(a2x_am2(lsize_am, nfields_sh_ar)) ! this will be obtained from atm instance
       if(iamroot_CPLID) then
          write(logunit,*) subname,' sharedFieldsAtmRof=', trim(sharedFieldsAtmRof)
           write(logunit,*) subname,' number of fields shared=', nfields_sh_ar
          write(logunit,*) subname,' seq_flds_a2x_fields_to_rof=', trim(seq_flds_a2x_fields_to_rof)
          write(logunit,*) subname,' lsize_am=', lsize_am
       endif
       a2racc_am(:,:) = 0.
       a2racc_am_cnt = 0

#endif
       allocate(a2r_rx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(a2r_rx(eri), rList=seq_flds_a2x_fields_to_rof, lsize=lsize_r)
          call mct_avect_zero(a2r_rx(eri))
       end do

       samegrid_ar = .true.
       if (trim(atm_gnam) /= trim(rof_gnam)) samegrid_ar = .false.

       if (atm_c2_rof) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fa2r'
          end if
          call seq_map_init_rcfile(mapper_Fa2r, atm(1), rof(1), &
               'seq_maps.rc','atm2rof_fmapname:','atm2rof_fmaptype:',samegrid_ar, &
               string='mapper_Fa2r initialization', esmf_map=esmf_map_flag)
! similar to a2o, prep_ocn 
#ifdef HAVE_MOAB
          ! Call moab intx only if atm  and river are init in moab
          if ((mbrxid .ge. 0) .and.  (mbaxid .ge. 0)) then
            appname = "ATM_ROF_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between rof and atm mesh
            idintx = 100*atm(1)%cplcompid + rof(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxar)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering atm rof intx'
              call shr_sys_abort(subname//' ERROR in registering atm rof intx')
            endif
            ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mbaxid, mbrxid, mbintxar)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in computing  atm rof intx'
              call shr_sys_abort(subname//' ERROR in computing atm rof intx')
            endif
            if (iamroot_CPLID) then
              write(logunit,*) 'iMOAB intersection between  atm and rof with id:', idintx
            end if
            ! we also need to compute the comm graph for the second hop, from the atm on coupler to the 
            ! atm for the intx atm-rof context (coverage)
            !    
            call seq_comm_getData(CPLID ,mpigrp=mpigrp_CPLID) 
            if (atm_pg_active) then
              type1 = 3; !  fv for both rof and atm; fv-cgll does not work anyway
            else
              type1 = 1 ! this does not work anyway in this direction
            endif
            type2 = 3;

            ierr = iMOAB_ComputeCommGraph( mbaxid, mbintxar, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        atm(1)%cplcompid, idintx)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing comm graph for second hop, rof-atm'
               call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, rof-atm')
            endif
            ! now take care of the mapper 
            if ( mapper_Fa2r%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fa2r%mbname) &
                             //' mapper_Fa2r'
                endif
            endif
            mapper_Fa2r%src_mbid = mbaxid
            mapper_Fa2r%tgt_mbid = mbrxid
            mapper_Fa2r%intx_mbid = mbintxar 
            mapper_Fa2r%src_context = atm(1)%cplcompid
            mapper_Fa2r%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fa2r%weight_identifier = wgtIdef
            mapper_Fa2r%mbname = 'mapper_Fa2r'
            ! because we will project fields from atm to rof grid, we need to define 
            ! rof a2x fields to rof grid on coupler side
            
            tagname = trim(seq_flds_a2x_fields_to_rof)//':norm8wt'//C_NULL_CHAR
            tagtype = 1 ! dense
            numco = 1 ! 
            ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags for seq_flds_a2x_fields_to_rof on rof cpl'
               call shr_sys_abort(subname//' ERROR in  defining tags for seq_flds_a2x_fields_to_rof on rof cpl')
            endif
            volumetric = 0 ! can be 1 only for FV->DGLL or FV->CGLL; 
            
            if (atm_pg_active) then
              dm1 = "fv"//C_NULL_CHAR
              dofnameS="GLOBAL_ID"//C_NULL_CHAR
              orderS = 1 !  fv-fv
            else ! this part does not work, anyway
              dm1 = "cgll"//C_NULL_CHAR
              dofnameS="GLOBAL_DOFS"//C_NULL_CHAR
              orderS = 4 ! np !  it should be 4
            endif
            dm2 = "fv"//C_NULL_CHAR
            dofnameT="GLOBAL_ID"//C_NULL_CHAR
            orderT = 1  !  not much arguing
            fNoBubble = 1
            monotonicity = 0 !
            noConserve = 0
            validate = 0 ! less verbose
            fInverseDistanceMap = 0
            if (iamroot_CPLID) then
               write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxar=', mbintxar, ' wgtIdef=', wgtIdef, &
                   'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameS), trim(dofnameT)
            endif
            ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxar, wgtIdef, &
                                               trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameS), trim(dofnameT) )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing ar weights '
               call shr_sys_abort(subname//' ERROR in computing ar weights ')
            endif

#ifdef MOABDEBUG
            wopts = C_NULL_CHAR
            call shr_mpi_commrank( mpicom_CPLID, rank )
            if (rank .lt. 5) then
              write(lnum,"(I0.2)")rank !
              outfile = 'intx_ar_'//trim(lnum)// '.h5m' // C_NULL_CHAR
              ierr = iMOAB_WriteMesh(mbintxar, outfile, wopts) ! write local intx file
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in writing intx ar file '
                call shr_sys_abort(subname//' ERROR in writing intx ar file ')
              endif
            endif
#endif
         end if ! if ((mbrxid .ge. 0) .and.  (mbaxid .ge. 0))
! endif HAVE_MOAB 
#endif

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2r'
          end if
          call seq_map_init_rcfile(mapper_Sa2r, atm(1), rof(1), &
               'seq_maps.rc','atm2rof_smapname:','atm2rof_smaptype:',samegrid_ar, &
               string='mapper_Sa2r initialization', esmf_map=esmf_map_flag)
#ifdef HAVE_MOAB
            ! now take care of the mapper, use the same one as before
            if ( mapper_Sa2r%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Sa2r%mbname) &
                             //' mapper_Sa2r'
                endif
            endif
            mapper_Sa2r%src_mbid = mbaxid
            mapper_Sa2r%tgt_mbid = mbrxid
            mapper_Sa2r%intx_mbid = mbintxar 
            mapper_Sa2r%src_context = atm(1)%cplcompid
            mapper_Sa2r%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Sa2r%weight_identifier = wgtIdef
            mapper_Sa2r%mbname = 'mapper_Sa2r'
#endif
       endif
	   
       call shr_sys_flush(logunit)

    end if

    if (rof_present .and. ocn_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       lsize_r = mct_aVect_lsize(x2r_rx)

       o2x_ox => component_get_c2x_cx(ocn(1))
       lsize_o = mct_aVect_lsize(o2x_ox)

       allocate(o2racc_ox(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_initSharedFields(o2x_ox, x2r_rx, o2racc_ox(eoi), lsize=lsize_o)
          call mct_aVect_zero(o2racc_ox(eoi))
       end do
       o2racc_ox_cnt = 0
#ifdef HAVE_MOAB

       ! this o2racc_om will be over ocn size 
       sharedFieldsOcnRof=''
       nfields_sh_or = mct_aVect_nRAttr(o2racc_ox(1))
       if ( nfields_sh_or /= 0 ) sharedFieldsOcnRof = trim( mct_aVect_exportRList2c(o2racc_ox(1)) )
       tagname = trim(sharedFieldsOcnRof)//C_NULL_CHAR
      
      ! find the size of ocn mesh locally
      ! find out the number of local elements in moab mesh ocn instance on coupler
      ierr  = iMOAB_GetMeshInfo ( mboxid, nvert, nvise, nbl, nsurf, nvisBC )
      if (ierr .ne. 0) then
         write(logunit,*) subname,' error in getting info '
         call shr_sys_abort(subname//' error in getting info ')
      endif
      ! ocn is fully cell now
      lsize_om = nvise(1)
      allocate(o2racc_om(lsize_om, nfields_sh_or))
      allocate(o2r_om2(lsize_om, nfields_sh_or)) ! this will be obtained from ocn instance
      if(iamroot_CPLID) then
         write(logunit,*) subname,' sharedFieldsOcnRof=', trim(sharedFieldsOcnRof)
         write(logunit,*) subname,' number of field shared ocn rof=',nfields_sh_or
         write(logunit,*) subname,' seq_flds_o2x_fields_to_rof=', trim(seq_flds_o2x_fields_to_rof)
         write(logunit,*) subname,' lsize_om=', lsize_om
      endif
      o2racc_om(:,:) = 0.
      o2racc_om_cnt = 0
#endif
       allocate(o2r_rx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(o2r_rx(eri), rList=seq_flds_o2x_fields_to_rof, lsize=lsize_r)
          call mct_avect_zero(o2r_rx(eri))
       end do

       samegrid_ro = .true.
       if (trim(ocn_gnam) /= trim(rof_gnam)) samegrid_ro = .false.

       if (ocn_c2_rof) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_So2r'
          end if
          call seq_map_init_rcfile(mapper_So2r, ocn(1), rof(1), &
               'seq_maps.rc','ocn2rof_smapname:','ocn2rof_smaptype:',samegrid_ro, &
               string='mapper_So2r initialization', esmf_map=esmf_map_flag)

       endif

       call shr_sys_flush(logunit)

    end if

  end subroutine prep_rof_init

  !================================================================================================
  subroutine prep_rof_accum_lnd(timer)

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
    character(*), parameter  :: subname = '(prep_rof_accum_lnd)'
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

  end subroutine prep_rof_accum_lnd

!================================================================================================
  subroutine prep_rof_accum_lnd_moab()

   use iMOAB , only :  iMOAB_GetDoubleTagStorage
   !---------------------------------------------------------------
   ! Description
   ! Accumulate land input to river component
   !
   !
   ! Local Variables
   character(CXX) ::tagname
   integer :: arrsize, ent_type, ierr
   character(*), parameter  :: subname = '(prep_rof_accum_lnd_moab)'
   !---------------------------------------------------------------

   ! do eli = 1,num_inst_lnd
   !    l2x_lx => component_get_c2x_cx(lnd(eli))
   !    if (l2racc_lx_cnt == 0) then
   !       call mct_avect_copy(l2x_lx, l2racc_lx(eli))
   !    else
   !       call mct_avect_accum(l2x_lx, l2racc_lx(eli))
   !    endif
   ! end do
   ! first, get l2x_lm2 from land coupler instance
   tagname = trim(sharedFieldsLndRof)//C_NULL_CHAR
   arrsize = nfields_sh_lr * lsize_lm
   ent_type = 1 ! cell type
   ierr = iMOAB_GetDoubleTagStorage ( mblxid, tagname, arrsize , ent_type, l2x_lm2)
   if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting shared fields from land instance ')
   endif
   ! big assumption is that l2x_lm2 is the same size as l2racc_lm
   if (l2racc_lm_cnt == 0) then
      l2racc_lm = l2x_lm2
   else
      l2racc_lm = l2racc_lm + l2x_lm2
   endif
   l2racc_lm_cnt = l2racc_lm_cnt + 1

 end subroutine prep_rof_accum_lnd_moab

  !================================================================================================

  subroutine prep_rof_accum_atm(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate atmosphere input to river component
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    type(mct_aVect), pointer :: a2x_ax
    character(*), parameter  :: subname = '(prep_rof_accum_atm)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)

    do eai = 1,num_inst_atm
       a2x_ax => component_get_c2x_cx(atm(eai))
       if (a2racc_ax_cnt == 0) then
          call mct_avect_copy(a2x_ax, a2racc_ax(eai))
       else
          call mct_avect_accum(a2x_ax, a2racc_ax(eai))
       endif
    end do
    a2racc_ax_cnt = a2racc_ax_cnt + 1

    call t_drvstopf (trim(timer))

  end subroutine prep_rof_accum_atm

!================================================================================================

  subroutine prep_rof_accum_atm_moab()

   use iMOAB , only :  iMOAB_GetDoubleTagStorage
   !---------------------------------------------------------------
   ! Description
   ! Accumulate atmosphere input to river component
   !
   !
   ! Local Variables
   character(CXX) ::tagname
   integer :: arrsize, ent_type, ierr
   character(*), parameter  :: subname = '(prep_rof_accum_atm_moab)'
   !---------------------------------------------------------------

   tagname = trim(sharedFieldsAtmRof)//C_NULL_CHAR
   arrsize = nfields_sh_ar * lsize_am
   ent_type = 1 ! cell type
   ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, a2x_am2)
   if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting shared fields from atm instance ')
   endif
   ! big assumption is that a2x_am2 is the same size as a2racc_am
   if (a2racc_am_cnt == 0) then
      a2racc_am = a2x_am2
   else
      a2racc_am = a2racc_am + a2x_am2
   endif
   a2racc_am_cnt = a2racc_am_cnt + 1
   !---------------------------------------------------------------

   ! call t_drvstartf (trim(timer),barrier=mpicom_CPLID)

   ! do eai = 1,num_inst_atm
   !    a2x_ax => component_get_c2x_cx(atm(eai))
   !    if (a2racc_ax_cnt == 0) then
   !       call mct_avect_copy(a2x_ax, a2racc_ax(eai))
   !    else
   !       call mct_avect_accum(a2x_ax, a2racc_ax(eai))
   !    endif
   ! end do
   ! a2racc_ax_cnt = a2racc_ax_cnt + 1


 end subroutine prep_rof_accum_atm_moab

  !================================================================================================
  subroutine prep_rof_accum_ocn(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate ocean input to river component
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi
    type(mct_aVect), pointer :: o2x_ox
    character(*), parameter  :: subname = '(prep_rof_accum_ocn)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)

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

  end subroutine prep_rof_accum_ocn

subroutine prep_rof_accum_ocn_moab()

    !---------------------------------------------------------------
    ! Description
    ! Accumulate ocean input to river component
    !
    !
    ! Local Variables

use iMOAB , only :  iMOAB_GetDoubleTagStorage
   !---------------------------------------------------------------
   ! Description
   ! Accumulate atmosphere input to river component
   !
   !
   ! Local Variables
   character(CXX) ::tagname
   integer :: arrsize, ent_type, ierr
    character(*), parameter  :: subname = '(prep_rof_accum_ocn_moab)'
   !---------------------------------------------------------------

   tagname = trim(sharedFieldsOcnRof)//C_NULL_CHAR
   arrsize = nfields_sh_or * lsize_om
   ent_type = 1 ! cell type
   ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, o2r_om2)
   if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting shared fields from ocn instance ')
   endif
   ! big assumption is that o2r_om2 is the same size as o2racc_om
   if (o2racc_om_cnt == 0) then
      o2racc_om = o2r_om2
   else
      o2racc_om = o2racc_om + o2r_om2
   endif
   o2racc_om_cnt = o2racc_om_cnt + 1


    !---------------------------------------------------------------


   !  do eoi = 1,num_inst_ocn
   !     o2x_ox => component_get_c2x_cx(ocn(eoi))
   !     if (o2racc_ox_cnt == 0) then
   !        call mct_avect_copy(o2x_ox, o2racc_ox(eoi))
   !     else
   !        call mct_avect_accum(o2x_ox, o2racc_ox(eoi))
   !     endif
   !  end do
   !  o2racc_ox_cnt = o2racc_ox_cnt + 1


  end subroutine prep_rof_accum_ocn_moab

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
    integer :: eri, eli, eai, eoi
    character(*), parameter :: subname = '(prep_rof_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if(l2racc_lx_cnt > 1) then
       do eri = 1,num_inst_rof
          eli = mod((eri-1),num_inst_lnd) + 1
          call mct_avect_avg(l2racc_lx(eli),l2racc_lx_cnt)
       enddo
    endif
    l2racc_lx_cnt = 0

    if((a2racc_ax_cnt > 1) .and. rof_heat) then
       do eri = 1,num_inst_rof
          eai = mod((eri-1),num_inst_atm) + 1
          call mct_avect_avg(a2racc_ax(eai),a2racc_ax_cnt)
       enddo
    endif
    a2racc_ax_cnt = 0

    if(o2racc_ox_cnt > 1) then
       do eri = 1,num_inst_rof
          eoi = mod((eri-1),num_inst_ocn) + 1
          call mct_avect_avg(o2racc_ox(eoi),o2racc_ox_cnt)
       enddo
    endif
    o2racc_ox_cnt = 0

    call t_drvstopf (trim(timer))

  end subroutine prep_rof_accum_avg

 !================================================================================================

  subroutine prep_rof_accum_avg_moab()

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of land, atm, ocn input to river component
    use iMOAB, only : iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
    use seq_comm_mct, only : num_moab_exports ! for debug
    ! Arguments
    !
    ! Local Variables
    character(CXX) ::tagname
    integer :: arrsize, ent_type, ierr
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
    character(*), parameter :: subname = '(prep_rof_accum_avg_moab)'
    !---------------------------------------------------------------
    if(l2racc_lm_cnt > 1) then
       l2racc_lm = 1./l2racc_lm_cnt*l2racc_lm
    endif
    l2racc_lm_cnt = 0
    ! set now the accumulated fields on land instance
    tagname = trim(sharedFieldsLndRof)//C_NULL_CHAR
    arrsize = nfields_sh_lr * lsize_lm
    ent_type = 1 ! cell type
    if (arrsize > 0) then
      ierr = iMOAB_SetDoubleTagStorage ( mblxid, tagname, arrsize , ent_type, l2racc_lm)
      if (ierr .ne. 0) then
         call shr_sys_abort(subname//' error in setting accumulated shared fields on rof on land instance ')
      endif
   endif

#ifdef MOABDEBUG
    if (mblxid .ge. 0 ) then !  we are on coupler pes, for sure
     write(lnum,"(I0.2)")num_moab_exports
     outfile = 'LndCplRofAvg'//trim(lnum)//'.h5m'//C_NULL_CHAR
     wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
     ierr = iMOAB_WriteMesh(mblxid, trim(outfile), trim(wopts))
     if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in writing land at rof accum ')
     endif
    endif
#endif

    if((a2racc_am_cnt > 1) .and. rof_heat) then
       a2racc_am = 1./a2racc_am_cnt * a2racc_am
    endif
    a2racc_am_cnt = 0
    ! set now the accumulated fields on atm instance
    tagname = trim(sharedFieldsAtmRof)//C_NULL_CHAR
    arrsize = nfields_sh_ar * lsize_am
    ent_type = 1 ! cell type
    if (arrsize > 0) then
      ierr = iMOAB_SetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, a2racc_am)
      if (ierr .ne. 0) then
         call shr_sys_abort(subname//' error in setting accumulated shared fields on rof on atm instance ')
      endif
   endif
#ifdef MOABDEBUG
    if (mbaxid .ge. 0 ) then !  we are on coupler pes, for sure
     write(lnum,"(I0.2)")num_moab_exports
     outfile = 'AtmCplRofAvg'//trim(lnum)//'.h5m'//C_NULL_CHAR
     wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
     ierr = iMOAB_WriteMesh(mbaxid, trim(outfile), trim(wopts))
     if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in writing atm at rof accum ')
     endif
    endif
#endif
    if(o2racc_om_cnt > 1) then
       o2racc_om = 1./o2racc_om_cnt *o2racc_om
    endif
    o2racc_om_cnt = 0
    ! set now the accumulated fields on ocn instance
    tagname = trim(sharedFieldsOcnRof)//C_NULL_CHAR
    arrsize = nfields_sh_or * lsize_om
    ent_type = 1 ! cell type
    if (arrsize > 0 ) then
      ierr = iMOAB_SetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, o2racc_om)
      if (ierr .ne. 0) then
         call shr_sys_abort(subname//' error in setting accumulated shared fields on rof on ocn instance ')
      endif
   endif
#ifdef MOABDEBUG
    if (mboxid .ge. 0 ) then !  we are on coupler pes, for sure
     write(lnum,"(I0.2)")num_moab_exports
     outfile = 'OcnCplRofAvg'//trim(lnum)//'.h5m'//C_NULL_CHAR
     wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
     ierr = iMOAB_WriteMesh(mboxid, trim(outfile), trim(wopts))
     if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in writing ocn at rof accum ')
     endif
    endif
#endif

  end subroutine prep_rof_accum_avg_moab


  !================================================================================================

  subroutine prep_rof_mrg(infodata, fractions_rx, timer_mrg, cime_model)

    !---------------------------------------------------------------
    ! Description
    ! Merge rof inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_rx(:)
    character(len=*)        , intent(in)    :: timer_mrg
    character(len=*)        , intent(in)    :: cime_model
    !
    ! Local Variables
    integer                  :: eri, efi, eoi
    type(mct_aVect), pointer :: x2r_rx
    character(*), parameter  :: subname = '(prep_rof_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg), barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       efi = mod((eri-1),num_inst_frc) + 1

       x2r_rx => component_get_x2c_cx(rof(eri))  ! This is actually modifying x2r_rx
       if(ocn_rof_two_way) then 
         call prep_rof_merge(l2r_rx(eri), a2r_rx(eri), fractions_rx(efi), x2r_rx, cime_model, o2x_r=o2r_rx(eri))
       else
         call prep_rof_merge(l2r_rx(eri), a2r_rx(eri), fractions_rx(efi), x2r_rx, cime_model)
       end if

    end do
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_rof_mrg

  !================================================================================================

  subroutine prep_rof_merge(l2x_r, a2x_r, fractions_r, x2r_r, cime_model,o2x_r)

    !-----------------------------------------------------------------------
    ! Description
    ! Merge land rof and ice forcing for rof input
    !
    ! Arguments
    type(mct_aVect),intent(in)    :: l2x_r
    type(mct_aVect),intent(in)    :: a2x_r
    type(mct_aVect),intent(in)    :: fractions_r
    type(mct_aVect),intent(inout) :: x2r_r
    character(len=*)        , intent(in)    :: cime_model
    type(mct_aVect),intent(in),optional  :: o2x_r
    !
    ! Local variables
    integer       :: i
    integer, save :: index_l2x_Flrl_rofsur
    integer, save :: index_l2x_Flrl_rofgwl
    integer, save :: index_l2x_Flrl_rofsub
    integer, save :: index_l2x_Flrl_rofdto
    integer, save :: index_l2x_Flrl_rofi
    integer, save :: index_l2x_Flrl_demand
    integer, save :: index_l2x_Flrl_irrig
    integer, save :: index_x2r_Flrl_rofsur
    integer, save :: index_x2r_Flrl_rofgwl
    integer, save :: index_x2r_Flrl_rofsub
    integer, save :: index_x2r_Flrl_rofdto
    integer, save :: index_x2r_Flrl_rofi
    integer, save :: index_x2r_Flrl_demand
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

    integer, save :: index_l2x_Flrl_Tqsur
    integer, save :: index_l2x_Flrl_Tqsub
    integer, save :: index_a2x_Sa_tbot
    integer, save :: index_a2x_Sa_pbot
    integer, save :: index_a2x_Sa_u
    integer, save :: index_a2x_Sa_v
    integer, save :: index_a2x_Sa_shum
    integer, save :: index_a2x_Faxa_swndr
    integer, save :: index_a2x_Faxa_swndf
    integer, save :: index_a2x_Faxa_swvdr
    integer, save :: index_a2x_Faxa_swvdf
    integer, save :: index_a2x_Faxa_lwdn
    integer, save :: index_x2r_Flrl_Tqsur
    integer, save :: index_x2r_Flrl_Tqsub
    integer, save :: index_x2r_Sa_tbot
    integer, save :: index_x2r_Sa_pbot
    integer, save :: index_x2r_Sa_u
    integer, save :: index_x2r_Sa_v
    integer, save :: index_x2r_Sa_shum
    integer, save :: index_x2r_Faxa_swndr
    integer, save :: index_x2r_Faxa_swndf
    integer, save :: index_x2r_Faxa_swvdr
    integer, save :: index_x2r_Faxa_swvdf
    integer, save :: index_x2r_Faxa_lwdn

    integer, save :: index_l2x_Flrl_inundinf
    integer, save :: index_x2r_Flrl_inundinf
    integer, save :: index_x2r_So_ssh
    integer, save :: index_o2x_So_ssh
    
    integer, save :: index_l2x_coszen_str
    integer, save :: index_x2r_coszen_str

    integer, save :: index_frac
    real(R8)      :: frac
    character(CL) :: fracstr
    logical, save :: first_time = .true.
    logical, save :: flds_wiso_rof = .false.
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
       if(trim(cime_model) .eq. 'e3sm') then
          index_l2x_Flrl_demand = mct_aVect_indexRA(l2x_r,'Flrl_demand' )
          index_x2r_Flrl_demand = mct_aVect_indexRA(x2r_r,'Flrl_demand' )
       endif
       index_x2r_Flrl_rofsur = mct_aVect_indexRA(x2r_r,'Flrl_rofsur' )
       index_x2r_Flrl_rofgwl = mct_aVect_indexRA(x2r_r,'Flrl_rofgwl' )
       index_x2r_Flrl_rofsub = mct_aVect_indexRA(x2r_r,'Flrl_rofsub' )
       index_x2r_Flrl_rofdto = mct_aVect_indexRA(x2r_r,'Flrl_rofdto' )
       index_x2r_Flrl_rofi   = mct_aVect_indexRA(x2r_r,'Flrl_rofi' )
       if (have_irrig_field) then
          index_x2r_Flrl_irrig  = mct_aVect_indexRA(x2r_r,'Flrl_irrig' )
       end if
       if(trim(cime_model) .eq. 'e3sm') then
         index_l2x_Flrl_Tqsur = mct_aVect_indexRA(l2x_r,'Flrl_Tqsur' )
         index_l2x_Flrl_Tqsub = mct_aVect_indexRA(l2x_r,'Flrl_Tqsub' )
         index_x2r_Flrl_Tqsur = mct_aVect_indexRA(x2r_r,'Flrl_Tqsur' )
         index_x2r_Flrl_Tqsub = mct_aVect_indexRA(x2r_r,'Flrl_Tqsub' )
       endif

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

       if (samegrid_al) then
          index_frac = mct_aVect_indexRA(fractions_r,"lfrac")
          fracstr = 'lfrac'
       else
          index_frac = mct_aVect_indexRA(fractions_r,"lfrin")
          fracstr = 'lfrin'
       endif

       mrgstr(index_x2r_Flrl_rofsur) = trim(mrgstr(index_x2r_Flrl_rofsur))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofsur'
       mrgstr(index_x2r_Flrl_rofgwl) = trim(mrgstr(index_x2r_Flrl_rofgwl))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofgwl'
       mrgstr(index_x2r_Flrl_rofsub) = trim(mrgstr(index_x2r_Flrl_rofsub))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofsub'
       mrgstr(index_x2r_Flrl_rofdto) = trim(mrgstr(index_x2r_Flrl_rofdto))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofdto'
       mrgstr(index_x2r_Flrl_rofi) = trim(mrgstr(index_x2r_Flrl_rofi))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofi'
       if (trim(cime_model).eq.'e3sm') then
          mrgstr(index_x2r_Flrl_demand) = trim(mrgstr(index_x2r_Flrl_demand))//' = '// &
               trim(fracstr)//'*l2x%Flrl_demand'
       endif
       if (have_irrig_field) then
          mrgstr(index_x2r_Flrl_irrig) = trim(mrgstr(index_x2r_Flrl_irrig))//' = '// &
               trim(fracstr)//'*l2x%Flrl_irrig'
       end if
       if(trim(cime_model) .eq. 'e3sm') then
          mrgstr(index_x2r_Flrl_Tqsur) = trim(mrgstr(index_x2r_Flrl_Tqsur))//' = '//'l2x%Flrl_Tqsur'
          mrgstr(index_x2r_Flrl_Tqsur) = trim(mrgstr(index_x2r_Flrl_Tqsub))//' = '//'l2x%Flrl_Tqsub'
       endif
       if ( flds_wiso_rof ) then
          mrgstr(index_x2r_Flrl_rofl_16O) = trim(mrgstr(index_x2r_Flrl_rofl_16O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofl_16O'
          mrgstr(index_x2r_Flrl_rofi_16O) = trim(mrgstr(index_x2r_Flrl_rofi_16O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofi_16O'
          mrgstr(index_x2r_Flrl_rofl_18O) = trim(mrgstr(index_x2r_Flrl_rofl_18O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofl_18O'
          mrgstr(index_x2r_Flrl_rofi_18O) = trim(mrgstr(index_x2r_Flrl_rofi_18O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofi_18O'
          mrgstr(index_x2r_Flrl_rofl_HDO) = trim(mrgstr(index_x2r_Flrl_rofl_HDO))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofl_HDO'
          mrgstr(index_x2r_Flrl_rofi_HDO) = trim(mrgstr(index_x2r_Flrl_rofi_HDO))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofi_HDO'
       end if
   
       if ( rof_heat ) then
          index_a2x_Sa_tbot    = mct_aVect_indexRA(a2x_r,'Sa_tbot')
          index_a2x_Sa_pbot    = mct_aVect_indexRA(a2x_r,'Sa_pbot')
          index_a2x_Sa_u       = mct_aVect_indexRA(a2x_r,'Sa_u')
          index_a2x_Sa_v       = mct_aVect_indexRA(a2x_r,'Sa_v')
          index_a2x_Sa_shum    = mct_aVect_indexRA(a2x_r,'Sa_shum')
          index_a2x_Faxa_swndr = mct_aVect_indexRA(a2x_r,'Faxa_swndr')
          index_a2x_Faxa_swndf = mct_aVect_indexRA(a2x_r,'Faxa_swndf')
          index_a2x_Faxa_swvdr = mct_aVect_indexRA(a2x_r,'Faxa_swvdr')
          index_a2x_Faxa_swvdf = mct_aVect_indexRA(a2x_r,'Faxa_swvdf')
          index_a2x_Faxa_lwdn  = mct_aVect_indexRA(a2x_r,'Faxa_lwdn')
     
          index_x2r_Sa_tbot    = mct_aVect_indexRA(x2r_r,'Sa_tbot')
          index_x2r_Sa_pbot    = mct_aVect_indexRA(x2r_r,'Sa_pbot')
          index_x2r_Sa_u       = mct_aVect_indexRA(x2r_r,'Sa_u')
          index_x2r_Sa_v       = mct_aVect_indexRA(x2r_r,'Sa_v')
          index_x2r_Sa_shum    = mct_aVect_indexRA(x2r_r,'Sa_shum')
          index_x2r_Faxa_swndr = mct_aVect_indexRA(x2r_r,'Faxa_swndr')
          index_x2r_Faxa_swndf = mct_aVect_indexRA(x2r_r,'Faxa_swndf')
          index_x2r_Faxa_swvdr = mct_aVect_indexRA(x2r_r,'Faxa_swvdr')
          index_x2r_Faxa_swvdf = mct_aVect_indexRA(x2r_r,'Faxa_swvdf')
          index_x2r_Faxa_lwdn  = mct_aVect_indexRA(x2r_r,'Faxa_lwdn')

          mrgstr(index_x2r_Sa_tbot)    = trim(mrgstr(index_x2r_Sa_tbot))//' = '//'a2x%Sa_tbot'
          mrgstr(index_x2r_Sa_pbot)    = trim(mrgstr(index_x2r_Sa_pbot))//' = '//'a2x%Sa_pbot'
          mrgstr(index_x2r_Sa_u)       = trim(mrgstr(index_x2r_Sa_u))//' = '//'a2x%Sa_u'
          mrgstr(index_x2r_Sa_v)       = trim(mrgstr(index_x2r_Sa_v))//' = '//'a2x%Sa_v'
          mrgstr(index_x2r_Sa_shum)    = trim(mrgstr(index_x2r_Sa_shum))//' = '//'a2x%Sa_shum'
          mrgstr(index_x2r_Faxa_swndr) = trim(mrgstr(index_x2r_Faxa_swndr))//' = '//'a2x%Faxa_swndr'
          mrgstr(index_x2r_Faxa_swndf) = trim(mrgstr(index_x2r_Faxa_swndf))//' = '//'a2x%Faxa_swndf'
          mrgstr(index_x2r_Faxa_swvdr) = trim(mrgstr(index_x2r_Faxa_swvdr))//' = '//'a2x%Faxa_swvdr'
          mrgstr(index_x2r_Faxa_swvdf) = trim(mrgstr(index_x2r_Faxa_swvdf))//' = '//'a2x%Faxa_swvdf'
          mrgstr(index_x2r_Faxa_lwdn)  = trim(mrgstr(index_x2r_Faxa_lwdn))//' = '//'a2x%Faxa_lwdn'

          if (lnd_rof_two_way) then
             index_l2x_Flrl_inundinf = mct_aVect_indexRA(l2x_r,'Flrl_inundinf')
             index_x2r_Flrl_inundinf = mct_aVect_indexRA(x2r_r,'Flrl_inundinf')
             mrgstr(index_x2r_Flrl_inundinf) = trim(mrgstr(index_x2r_Flrl_inundinf))//' = '//'l2x%Flrl_inundinf'
          endif

       endif 

       if (ocn_rof_two_way) then
          index_o2x_So_ssh = mct_aVect_indexRA(o2x_r,'So_ssh')
          index_x2r_So_ssh = mct_aVect_indexRA(x2r_r,'So_ssh')
          mrgstr(index_x2r_So_ssh)       = trim(mrgstr(index_x2r_So_ssh))//' = '//'o2x%So_ssh'
       endif

    end if

    do i = 1,lsize
       frac = fractions_r%rAttr(index_frac,i)
       x2r_r%rAttr(index_x2r_Flrl_rofsur,i) = l2x_r%rAttr(index_l2x_Flrl_rofsur,i) * frac
       x2r_r%rAttr(index_x2r_Flrl_rofgwl,i) = l2x_r%rAttr(index_l2x_Flrl_rofgwl,i) * frac
       x2r_r%rAttr(index_x2r_Flrl_rofsub,i) = l2x_r%rAttr(index_l2x_Flrl_rofsub,i) * frac
       x2r_r%rAttr(index_x2r_Flrl_rofdto,i) = l2x_r%rAttr(index_l2x_Flrl_rofdto,i) * frac
       x2r_r%rAttr(index_x2r_Flrl_rofi,i) = l2x_r%rAttr(index_l2x_Flrl_rofi,i) * frac
       if (trim(cime_model).eq.'e3sm') then
          x2r_r%rAttr(index_x2r_Flrl_demand,i) = l2x_r%rAttr(index_l2x_Flrl_demand,i) * frac
       endif
       if (have_irrig_field) then
          x2r_r%rAttr(index_x2r_Flrl_irrig,i) = l2x_r%rAttr(index_l2x_Flrl_irrig,i) * frac
       end if
       if(trim(cime_model) .eq. 'e3sm') then
         x2r_r%rAttr(index_x2r_Flrl_Tqsur,i) = l2x_r%rAttr(index_l2x_Flrl_Tqsur,i)
         x2r_r%rAttr(index_x2r_Flrl_Tqsub,i) = l2x_r%rAttr(index_l2x_Flrl_Tqsub,i)
       endif
       if ( flds_wiso_rof ) then
          x2r_r%rAttr(index_x2r_Flrl_rofl_16O,i) = l2x_r%rAttr(index_l2x_Flrl_rofl_16O,i) * frac
          x2r_r%rAttr(index_x2r_Flrl_rofi_16O,i) = l2x_r%rAttr(index_l2x_Flrl_rofi_16O,i) * frac
          x2r_r%rAttr(index_x2r_Flrl_rofl_18O,i) = l2x_r%rAttr(index_l2x_Flrl_rofl_18O,i) * frac
          x2r_r%rAttr(index_x2r_Flrl_rofi_18O,i) = l2x_r%rAttr(index_l2x_Flrl_rofi_18O,i) * frac
          x2r_r%rAttr(index_x2r_Flrl_rofl_HDO,i) = l2x_r%rAttr(index_l2x_Flrl_rofl_HDO,i) * frac
          x2r_r%rAttr(index_x2r_Flrl_rofi_HDO,i) = l2x_r%rAttr(index_l2x_Flrl_rofi_HDO,i) * frac
       end if
      
       if ( rof_heat ) then
          x2r_r%rAttr(index_x2r_Sa_tbot,i)    = a2x_r%rAttr(index_a2x_Sa_tbot,i)
          x2r_r%rAttr(index_x2r_Sa_pbot,i)    = a2x_r%rAttr(index_a2x_Sa_pbot,i)
          x2r_r%rAttr(index_x2r_Sa_u,i)       = a2x_r%rAttr(index_a2x_Sa_u,i)
          x2r_r%rAttr(index_x2r_Sa_v,i)       = a2x_r%rAttr(index_a2x_Sa_v,i)
          x2r_r%rAttr(index_x2r_Sa_shum,i)    = a2x_r%rAttr(index_a2x_Sa_shum,i)
          x2r_r%rAttr(index_x2r_Faxa_swndr,i) = a2x_r%rAttr(index_a2x_Faxa_swndr,i)
          x2r_r%rAttr(index_x2r_Faxa_swndf,i) = a2x_r%rAttr(index_a2x_Faxa_swndf,i)
          x2r_r%rAttr(index_x2r_Faxa_swvdr,i) = a2x_r%rAttr(index_a2x_Faxa_swvdr,i)
          x2r_r%rAttr(index_x2r_Faxa_swvdf,i) = a2x_r%rAttr(index_a2x_Faxa_swvdf,i)
          x2r_r%rAttr(index_x2r_Faxa_lwdn,i)  = a2x_r%rAttr(index_a2x_Faxa_lwdn,i)
       endif

       if (lnd_rof_two_way) then
         x2r_r%rAttr(index_x2r_Flrl_inundinf,i) = l2x_r%rAttr(index_l2x_Flrl_inundinf,i)
       endif

    end do

    if (ocn_rof_two_way) then
       do i =1,lsize
          x2r_r%rAttr(index_x2r_So_ssh,i)        = o2x_r%rAttr(index_o2x_So_ssh,i)
       enddo
    endif

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
#ifdef HAVE_MOAB
  subroutine prep_rof_mrg_moab  (infodata, cime_model)
   use iMOAB , only : iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, &
     iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
     use seq_comm_mct, only : num_moab_exports ! for debug

    type(seq_infodata_type) , intent(in)    :: infodata

    ! type(mct_aVect)         , intent(in)    :: fractions_rx(:) they should have been saved as tags on rof coupler component
    character(len=*)        , intent(in)    :: cime_model
    !-----------------------------------------------------------------------
    ! Description
    ! Merge land rof and ice forcing for rof input
    !
    ! used for indexing 
    type(mct_avect) , pointer   :: l2x_r
    type(mct_avect) , pointer   :: a2x_r
    type(mct_avect) , pointer    :: fractions_r
    type(mct_avect) , pointer :: x2r_r
    !
    ! Local variables
    integer       :: i
    integer, save :: index_l2x_Flrl_rofsur
    integer, save :: index_l2x_Flrl_rofgwl
    integer, save :: index_l2x_Flrl_rofsub
    integer, save :: index_l2x_Flrl_rofdto
    integer, save :: index_l2x_Flrl_rofi
    integer, save :: index_l2x_Flrl_demand
    integer, save :: index_l2x_Flrl_irrig
    integer, save :: index_x2r_Flrl_rofsur
    integer, save :: index_x2r_Flrl_rofgwl
    integer, save :: index_x2r_Flrl_rofsub
    integer, save :: index_x2r_Flrl_rofdto
    integer, save :: index_x2r_Flrl_rofi
    integer, save :: index_x2r_Flrl_demand
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
	
    integer, save :: index_l2x_Flrl_Tqsur
    integer, save :: index_l2x_Flrl_Tqsub
    integer, save :: index_a2x_Sa_tbot
    integer, save :: index_a2x_Sa_pbot
    integer, save :: index_a2x_Sa_u
    integer, save :: index_a2x_Sa_v
    integer, save :: index_a2x_Sa_shum
    integer, save :: index_a2x_Faxa_swndr
    integer, save :: index_a2x_Faxa_swndf
    integer, save :: index_a2x_Faxa_swvdr
    integer, save :: index_a2x_Faxa_swvdf
    integer, save :: index_a2x_Faxa_lwdn
    integer, save :: index_x2r_Flrl_Tqsur
    integer, save :: index_x2r_Flrl_Tqsub
    integer, save :: index_x2r_Sa_tbot
    integer, save :: index_x2r_Sa_pbot
    integer, save :: index_x2r_Sa_u
    integer, save :: index_x2r_Sa_v
    integer, save :: index_x2r_Sa_shum
    integer, save :: index_x2r_Faxa_swndr
    integer, save :: index_x2r_Faxa_swndf
    integer, save :: index_x2r_Faxa_swvdr
    integer, save :: index_x2r_Faxa_swvdf
    integer, save :: index_x2r_Faxa_lwdn
    
    integer, save :: index_l2x_coszen_str
    integer, save :: index_x2r_coszen_str

    integer, save :: index_frac
    real(R8)      :: frac
    character(CL) :: fracstr
    logical, save :: first_time = .true.
    logical, save :: flds_wiso_rof = .false.
    integer, save :: nflds, lsize
    logical       :: iamroot
    character(CL) :: field        ! field string
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(*), parameter   :: subname = '(prep_rof_mrg_moab) '

 integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
   
    character(CXX) ::tagname, mct_field
    integer :: ent_type, ierr, arrsize
    integer, save :: naflds, nlflds ! these are saved the first time 
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
#ifdef MOABCOMP
    real(R8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
#endif 
    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    
! character(*),parameter :: fraclist_r = 'lfrac:lfrin:rfrac' 
    if (first_time) then
      ! find out the number of local elements in moab mesh rof instance on coupler
      ierr  = iMOAB_GetMeshInfo ( mbrxid, nvert, nvise, nbl, nsurf, nvisBC )
      if (ierr .ne. 0) then
            write(logunit,*) subname,' error in getting info '
            call shr_sys_abort(subname//' error in getting info ')
      endif
      lsize = nvise(1) ! number of active cells
       ! mct avs are used just for their fields metadata, not the actual reals 
       ! (name of the fields)
       ! need these always, not only the first time
      l2x_r => l2r_rx(1)
      a2x_r => a2r_rx(1)

      x2r_r => component_get_x2c_cx(rof(1))
      nflds = mct_aVect_nRattr(x2r_r) ! these are saved after first time
      naflds = mct_aVect_nRattr(a2x_r)
      nlflds = mct_aVect_nRattr(l2x_r)

      allocate(x2r_rm (lsize, nflds))
      allocate(a2x_rm (lsize, naflds))
      allocate(l2x_rm (lsize, nlflds))
      ! allocate fractions too
      ! use the fraclist fraclist_r = 'lfrac:lfrin:rfrac'
      allocate(fractions_rm(lsize,3)) ! there are 3 fields here

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
       if(trim(cime_model) .eq. 'e3sm') then
          index_l2x_Flrl_demand = mct_aVect_indexRA(l2x_r,'Flrl_demand' )
          index_x2r_Flrl_demand = mct_aVect_indexRA(x2r_r,'Flrl_demand' )
       endif
       index_x2r_Flrl_rofsur = mct_aVect_indexRA(x2r_r,'Flrl_rofsur' )
       index_x2r_Flrl_rofgwl = mct_aVect_indexRA(x2r_r,'Flrl_rofgwl' )
       index_x2r_Flrl_rofsub = mct_aVect_indexRA(x2r_r,'Flrl_rofsub' )
       index_x2r_Flrl_rofdto = mct_aVect_indexRA(x2r_r,'Flrl_rofdto' )
       index_x2r_Flrl_rofi   = mct_aVect_indexRA(x2r_r,'Flrl_rofi' )
       if (have_irrig_field) then
          index_x2r_Flrl_irrig  = mct_aVect_indexRA(x2r_r,'Flrl_irrig' )
       end if
       if(trim(cime_model) .eq. 'e3sm') then
         index_l2x_Flrl_Tqsur = mct_aVect_indexRA(l2x_r,'Flrl_Tqsur' )
         index_l2x_Flrl_Tqsub = mct_aVect_indexRA(l2x_r,'Flrl_Tqsub' )
         index_x2r_Flrl_Tqsur = mct_aVect_indexRA(x2r_r,'Flrl_Tqsur' )
         index_x2r_Flrl_Tqsub = mct_aVect_indexRA(x2r_r,'Flrl_Tqsub' )
       endif

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

       if (samegrid_al) then ! fraclist_r = 'lfrac:lfrin:rfrac'
       ! check, in our case, is it always 1  ? 
          index_frac = 1 ! mct_aVect_indexRA(fractions_r,"lfrac")
          fracstr = 'lfrac'
       else
          index_frac = 2 ! mct_aVect_indexRA(fractions_r,"lfrin")
          fracstr = 'lfrin'
       endif

       mrgstr(index_x2r_Flrl_rofsur) = trim(mrgstr(index_x2r_Flrl_rofsur))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofsur'
       mrgstr(index_x2r_Flrl_rofgwl) = trim(mrgstr(index_x2r_Flrl_rofgwl))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofgwl'
       mrgstr(index_x2r_Flrl_rofsub) = trim(mrgstr(index_x2r_Flrl_rofsub))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofsub'
       mrgstr(index_x2r_Flrl_rofdto) = trim(mrgstr(index_x2r_Flrl_rofdto))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofdto'
       mrgstr(index_x2r_Flrl_rofi) = trim(mrgstr(index_x2r_Flrl_rofi))//' = '// &
            trim(fracstr)//'*l2x%Flrl_rofi'
       if (trim(cime_model).eq.'e3sm') then
          mrgstr(index_x2r_Flrl_demand) = trim(mrgstr(index_x2r_Flrl_demand))//' = '// &
               trim(fracstr)//'*l2x%Flrl_demand'
       endif
       if (have_irrig_field) then
          mrgstr(index_x2r_Flrl_irrig) = trim(mrgstr(index_x2r_Flrl_irrig))//' = '// &
               trim(fracstr)//'*l2x%Flrl_irrig'
       end if
       if(trim(cime_model) .eq. 'e3sm') then
          mrgstr(index_x2r_Flrl_Tqsur) = trim(mrgstr(index_x2r_Flrl_Tqsur))//' = '//'l2x%Flrl_Tqsur'
          mrgstr(index_x2r_Flrl_Tqsur) = trim(mrgstr(index_x2r_Flrl_Tqsub))//' = '//'l2x%Flrl_Tqsub'
       endif
       if ( flds_wiso_rof ) then
          mrgstr(index_x2r_Flrl_rofl_16O) = trim(mrgstr(index_x2r_Flrl_rofl_16O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofl_16O'
          mrgstr(index_x2r_Flrl_rofi_16O) = trim(mrgstr(index_x2r_Flrl_rofi_16O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofi_16O'
          mrgstr(index_x2r_Flrl_rofl_18O) = trim(mrgstr(index_x2r_Flrl_rofl_18O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofl_18O'
          mrgstr(index_x2r_Flrl_rofi_18O) = trim(mrgstr(index_x2r_Flrl_rofi_18O))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofi_18O'
          mrgstr(index_x2r_Flrl_rofl_HDO) = trim(mrgstr(index_x2r_Flrl_rofl_HDO))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofl_HDO'
          mrgstr(index_x2r_Flrl_rofi_HDO) = trim(mrgstr(index_x2r_Flrl_rofi_HDO))//' = '// &
               trim(fracstr)//'*l2x%Flrl_rofi_HDO'
       end if
	   
       if ( rof_heat ) then
          index_a2x_Sa_tbot    = mct_aVect_indexRA(a2x_r,'Sa_tbot')
          index_a2x_Sa_pbot    = mct_aVect_indexRA(a2x_r,'Sa_pbot')
          index_a2x_Sa_u       = mct_aVect_indexRA(a2x_r,'Sa_u')
          index_a2x_Sa_v       = mct_aVect_indexRA(a2x_r,'Sa_v')
          index_a2x_Sa_shum    = mct_aVect_indexRA(a2x_r,'Sa_shum')
          index_a2x_Faxa_swndr = mct_aVect_indexRA(a2x_r,'Faxa_swndr')
          index_a2x_Faxa_swndf = mct_aVect_indexRA(a2x_r,'Faxa_swndf')
          index_a2x_Faxa_swvdr = mct_aVect_indexRA(a2x_r,'Faxa_swvdr')
          index_a2x_Faxa_swvdf = mct_aVect_indexRA(a2x_r,'Faxa_swvdf')
          index_a2x_Faxa_lwdn  = mct_aVect_indexRA(a2x_r,'Faxa_lwdn')
     
          index_x2r_Sa_tbot    = mct_aVect_indexRA(x2r_r,'Sa_tbot')
          index_x2r_Sa_pbot    = mct_aVect_indexRA(x2r_r,'Sa_pbot')
          index_x2r_Sa_u       = mct_aVect_indexRA(x2r_r,'Sa_u')
          index_x2r_Sa_v       = mct_aVect_indexRA(x2r_r,'Sa_v')
          index_x2r_Sa_shum    = mct_aVect_indexRA(x2r_r,'Sa_shum')
          index_x2r_Faxa_swndr = mct_aVect_indexRA(x2r_r,'Faxa_swndr')
          index_x2r_Faxa_swndf = mct_aVect_indexRA(x2r_r,'Faxa_swndf')
          index_x2r_Faxa_swvdr = mct_aVect_indexRA(x2r_r,'Faxa_swvdr')
          index_x2r_Faxa_swvdf = mct_aVect_indexRA(x2r_r,'Faxa_swvdf')
          index_x2r_Faxa_lwdn  = mct_aVect_indexRA(x2r_r,'Faxa_lwdn')

          mrgstr(index_x2r_Sa_tbot)    = trim(mrgstr(index_x2r_Sa_tbot))//' = '//'a2x%Sa_tbot'
          mrgstr(index_x2r_Sa_pbot)    = trim(mrgstr(index_x2r_Sa_pbot))//' = '//'a2x%Sa_pbot'
          mrgstr(index_x2r_Sa_u)       = trim(mrgstr(index_x2r_Sa_u))//' = '//'a2x%Sa_u'
          mrgstr(index_x2r_Sa_v)       = trim(mrgstr(index_x2r_Sa_v))//' = '//'a2x%Sa_v'
          mrgstr(index_x2r_Sa_shum)    = trim(mrgstr(index_x2r_Sa_shum))//' = '//'a2x%Sa_shum'
          mrgstr(index_x2r_Faxa_swndr) = trim(mrgstr(index_x2r_Faxa_swndr))//' = '//'a2x%Faxa_swndr'
          mrgstr(index_x2r_Faxa_swndf) = trim(mrgstr(index_x2r_Faxa_swndf))//' = '//'a2x%Faxa_swndf'
          mrgstr(index_x2r_Faxa_swvdr) = trim(mrgstr(index_x2r_Faxa_swvdr))//' = '//'a2x%Faxa_swvdr'
          mrgstr(index_x2r_Faxa_swvdf) = trim(mrgstr(index_x2r_Faxa_swvdf))//' = '//'a2x%Faxa_swvdf'
          mrgstr(index_x2r_Faxa_lwdn)  = trim(mrgstr(index_x2r_Faxa_lwdn))//' = '//'a2x%Faxa_lwdn'
       endif 

    end if
! fill in with data from moab tags

! fill with fractions from river instance
! fractions_rm, 
    ent_type = 1 ! cells
    tagname = 'lfrac:lfrin:rfrac'//C_NULL_CHAR
    arrsize = 3 * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mbrxid, tagname, arrsize , ent_type, fractions_rm)
    if (ierr .ne. 0) then
         call shr_sys_abort(subname//' error in getting fractions_om from rof instance ')
    endif
   ! fill the r2x_rm, etc double array fields nflds
    tagname = trim(seq_flds_x2r_fields)//C_NULL_CHAR
    arrsize = nflds * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mbrxid, tagname, arrsize , ent_type, x2r_rm)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting x2r_rm array ')
    endif    
    ! a2x_rm (lsize, naflds))

    tagname = trim(seq_flds_a2x_fields_to_rof)//C_NULL_CHAR
    arrsize = naflds * lsize !        allocate (a2x_rm (lsize, naflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbrxid, tagname, arrsize , ent_type, a2x_rm)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting a2x_rm array ')
    endif
    !  l2x_rm
    tagname = trim(seq_flds_l2x_fluxes_to_rof)//C_NULL_CHAR
    arrsize = nlflds * lsize !        
    ierr = iMOAB_GetDoubleTagStorage ( mbrxid, tagname, arrsize , ent_type, l2x_rm)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting l2x_rm array ')
    endif

! replace x2r_r%rAttr(index ,i) with x2r_rm(i,index), etc from formula in prep_rof_merge
! x2r_r%rAttr( -> x2r_rm(i, %rAttr( -> m(,i)  ,i) -> )
    do i = 1,lsize
       frac = fractions_rm(i,index_frac)
       x2r_rm(i,index_x2r_Flrl_rofsur) = l2x_rm(i,index_l2x_Flrl_rofsur) * frac
       x2r_rm(i,index_x2r_Flrl_rofgwl) = l2x_rm(i,index_l2x_Flrl_rofgwl) * frac
       x2r_rm(i,index_x2r_Flrl_rofsub) = l2x_rm(i,index_l2x_Flrl_rofsub) * frac
       x2r_rm(i,index_x2r_Flrl_rofdto) = l2x_rm(i,index_l2x_Flrl_rofdto) * frac
       x2r_rm(i,index_x2r_Flrl_rofi) = l2x_rm(i,index_l2x_Flrl_rofi) * frac
       if (trim(cime_model).eq.'e3sm') then
          x2r_rm(i,index_x2r_Flrl_demand) = l2x_rm(i,index_l2x_Flrl_demand) * frac
       endif
       if (have_irrig_field) then
          x2r_rm(i,index_x2r_Flrl_irrig) = l2x_rm(i,index_l2x_Flrl_irrig) * frac
       end if
       if(trim(cime_model) .eq. 'e3sm') then
         x2r_rm(i,index_x2r_Flrl_Tqsur) = l2x_rm(i,index_l2x_Flrl_Tqsur)
         x2r_rm(i,index_x2r_Flrl_Tqsub) = l2x_rm(i,index_l2x_Flrl_Tqsub)
       endif
       if ( flds_wiso_rof ) then
          x2r_rm(i,index_x2r_Flrl_rofl_16O) = l2x_rm(i,index_l2x_Flrl_rofl_16O) * frac
          x2r_rm(i,index_x2r_Flrl_rofi_16O) = l2x_rm(i,index_l2x_Flrl_rofi_16O) * frac
          x2r_rm(i,index_x2r_Flrl_rofl_18O) = l2x_rm(i,index_l2x_Flrl_rofl_18O) * frac
          x2r_rm(i,index_x2r_Flrl_rofi_18O) = l2x_rm(i,index_l2x_Flrl_rofi_18O) * frac
          x2r_rm(i,index_x2r_Flrl_rofl_HDO) = l2x_rm(i,index_l2x_Flrl_rofl_HDO) * frac
          x2r_rm(i,index_x2r_Flrl_rofi_HDO) = l2x_rm(i,index_l2x_Flrl_rofi_HDO) * frac
       end if
      
       if ( rof_heat ) then
          x2r_rm(i,index_x2r_Sa_tbot)    = a2x_rm(i,index_a2x_Sa_tbot)
          x2r_rm(i,index_x2r_Sa_pbot)    = a2x_rm(i,index_a2x_Sa_pbot)
          x2r_rm(i,index_x2r_Sa_u)       = a2x_rm(i,index_a2x_Sa_u)
          x2r_rm(i,index_x2r_Sa_v)       = a2x_rm(i,index_a2x_Sa_v)
          x2r_rm(i,index_x2r_Sa_shum)    = a2x_rm(i,index_a2x_Sa_shum)
          x2r_rm(i,index_x2r_Faxa_swndr) = a2x_rm(i,index_a2x_Faxa_swndr)
          x2r_rm(i,index_x2r_Faxa_swndf) = a2x_rm(i,index_a2x_Faxa_swndf)
          x2r_rm(i,index_x2r_Faxa_swvdr) = a2x_rm(i,index_a2x_Faxa_swvdr)
          x2r_rm(i,index_x2r_Faxa_swvdf) = a2x_rm(i,index_a2x_Faxa_swvdf)
          x2r_rm(i,index_x2r_Faxa_lwdn)  = a2x_rm(i,index_a2x_Faxa_lwdn)
       endif

    end do
    ! after we are done, set x2r_rm to the mbrxid

    tagname = trim(seq_flds_x2r_fields)//C_NULL_CHAR
    arrsize = nflds * lsize
    ierr = iMOAB_SetDoubleTagStorage ( mbrxid, tagname, arrsize , ent_type, x2r_rm)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in setting x2r_rm array ')
    endif
    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,nflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif


#ifdef MOABCOMP
  !compare_mct_av_moab_tag(comp, attrVect, field, imoabApp, tag_name, ent_type, difference)
    x2r_r => component_get_x2c_cx(rof(1))
    ! loop over all fields in seq_flds_x2r_fields
    call mct_list_init(temp_list ,seq_flds_x2r_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! cell for river
    if (iamroot) print *, num_moab_exports, trim(seq_flds_x2r_fields)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call compare_mct_av_moab_tag(rof(1), x2r_r, mct_field,  mbrxid, tagname, ent_type, difference, first_time)
    enddo
    call mct_list_clean(temp_list)
#endif
    first_time = .false.
#ifdef MOABDEBUG

    if (mbrxid .ge. 0 ) then !  we are on coupler pes, for sure
     write(lnum,"(I0.2)")num_moab_exports
     outfile = 'RofCplAftMm'//trim(lnum)//'.h5m'//C_NULL_CHAR
     wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
     ierr = iMOAB_WriteMesh(mbrxid, trim(outfile), trim(wopts))
   endif
#endif

  end subroutine prep_rof_mrg_moab
#endif 
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

  subroutine prep_rof_calc_a2r_rx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create a2r_rx (note that a2r_rx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri, eai
    type(mct_avect), pointer :: r2x_rx
    character(*), parameter :: subname = '(prep_rof_calc_a2r_rx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eai = mod((eri-1),num_inst_atm) + 1
       r2x_rx => component_get_c2x_cx(rof(eri))
       call seq_map_map(mapper_Sa2r, a2racc_ax(eai), a2r_rx(eri), fldlist=seq_flds_a2x_states_to_rof, norm=.true.)
       call seq_map_map(mapper_Fa2r, a2racc_ax(eai), a2r_rx(eri), fldlist=seq_flds_a2x_fluxes_to_rof, norm=.true.)
    end do
    call t_drvstopf  (trim(timer))

  end subroutine prep_rof_calc_a2r_rx

  !================================================================================================
  subroutine prep_rof_calc_o2r_rx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2r_rx (note that o2r_rx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri, eoi
    type(mct_avect), pointer :: r2x_rx
    character(*), parameter :: subname = '(prep_rof_calc_o2r_rx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       eoi = mod((eri-1),num_inst_ocn) + 1
       r2x_rx => component_get_c2x_cx(rof(eri))
       call seq_map_map(mapper_So2r, o2racc_ox(eoi), o2r_rx(eri), fldlist=seq_flds_o2x_states_to_rof, norm=.true.)
    end do
    call t_drvstopf  (trim(timer))

  end subroutine prep_rof_calc_o2r_rx

  !================================================================================================

  function prep_rof_get_l2racc_lx()
    type(mct_aVect), pointer :: prep_rof_get_l2racc_lx(:)
    prep_rof_get_l2racc_lx => l2racc_lx(:)
  end function prep_rof_get_l2racc_lx

  function prep_rof_get_l2racc_lx_cnt()
    integer, pointer :: prep_rof_get_l2racc_lx_cnt
    prep_rof_get_l2racc_lx_cnt => l2racc_lx_cnt
  end function prep_rof_get_l2racc_lx_cnt

  function prep_rof_get_o2racc_ox()
    type(mct_aVect), pointer :: prep_rof_get_o2racc_ox(:)
    prep_rof_get_o2racc_ox => o2racc_ox(:)
  end function prep_rof_get_o2racc_ox

  function prep_rof_get_o2racc_ox_cnt()
    integer, pointer :: prep_rof_get_o2racc_ox_cnt
    prep_rof_get_o2racc_ox_cnt => o2racc_ox_cnt
  end function prep_rof_get_o2racc_ox_cnt

  function prep_rof_get_mapper_Fl2r()
    type(seq_map), pointer :: prep_rof_get_mapper_Fl2r
    prep_rof_get_mapper_Fl2r => mapper_Fl2r
  end function prep_rof_get_mapper_Fl2r

  function prep_rof_get_a2racc_ax()
    type(mct_aVect), pointer :: prep_rof_get_a2racc_ax(:)
    prep_rof_get_a2racc_ax => a2racc_ax(:)
  end function prep_rof_get_a2racc_ax

  function prep_rof_get_a2racc_ax_cnt()
    integer, pointer :: prep_rof_get_a2racc_ax_cnt
    prep_rof_get_a2racc_ax_cnt => a2racc_ax_cnt
  end function prep_rof_get_a2racc_ax_cnt

  function prep_rof_get_mapper_Fa2r()
    type(seq_map), pointer :: prep_rof_get_mapper_Fa2r
    prep_rof_get_mapper_Fa2r => mapper_Fa2r
  end function prep_rof_get_mapper_Fa2r

  function prep_rof_get_mapper_Sa2r()
    type(seq_map), pointer :: prep_rof_get_mapper_Sa2r
    prep_rof_get_mapper_Sa2r => mapper_Sa2r
  end function prep_rof_get_mapper_Sa2r

 ! moab 
 ! for ocean
  function prep_rof_get_o2racc_om_cnt()
    integer, pointer :: prep_rof_get_o2racc_om_cnt
    prep_rof_get_o2racc_om_cnt => o2racc_om_cnt
  end function prep_rof_get_o2racc_om_cnt

  function prep_rof_get_o2racc_om()
   real(R8), DIMENSION(:, :), pointer :: prep_rof_get_o2racc_om
   prep_rof_get_o2racc_om => o2racc_om
  end function prep_rof_get_o2racc_om

  function prep_rof_get_sharedFieldsOcnRof()
     character(CXX) :: prep_rof_get_sharedFieldsOcnRof
     prep_rof_get_sharedFieldsOcnRof = sharedFieldsOcnRof
  end function prep_rof_get_sharedFieldsOcnRof
 
  ! for land
  function prep_rof_get_l2racc_lm_cnt()
    integer, pointer :: prep_rof_get_l2racc_lm_cnt
    prep_rof_get_l2racc_lm_cnt => l2racc_lm_cnt
  end function prep_rof_get_l2racc_lm_cnt

  function prep_rof_get_l2racc_lm()
    real(R8), DIMENSION(:, :), pointer :: prep_rof_get_l2racc_lm
    prep_rof_get_l2racc_lm => l2racc_lm
  end function prep_rof_get_l2racc_lm

  function prep_rof_get_sharedFieldsLndRof()
     character(CXX) :: prep_rof_get_sharedFieldsLndRof
     prep_rof_get_sharedFieldsLndRof = sharedFieldsLndRof
  end function prep_rof_get_sharedFieldsLndRof


end module prep_rof_mod
