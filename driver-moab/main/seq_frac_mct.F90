! !MODULE: seq_frac_mct -- handles surface fractions.
!
!  Fraction Notes: tcraig, august 2008
!  Assumes is running on CPLID pes
!
!  the fractions fields are now afrac, ifrac, ofrac, lfrac, and lfrin.
!    afrac = fraction of atm on a grid
!    lfrac = fraction of lnd on a grid
!    ifrac = fraction of ice on a grid
!    ofrac = fraction of ocn on a grid
!    lfrin = land fraction defined by the land model
!    ifrad = fraction of ocn on a grid at last radiation time
!    ofrad = fraction of ice on a grid at last radiation time
!      afrac, lfrac, ifrac, and ofrac are the self-consistent values in the
!      system.  lfrin is the fraction on the land grid and is allowed to
!      vary from the self-consistent value as descibed below.  ifrad
!      and ofrad are needed for the swnet calculation.
!  the fractions fields are defined for each grid in the fraction bundles as
!    needed as follows.
!    character(*),parameter :: fraclist_a = 'afrac:ifrac:ofrac:lfrac:lfrin'
!    character(*),parameter :: fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
!    character(*),parameter :: fraclist_i = 'afrac:ifrac:ofrac'
!    character(*),parameter :: fraclist_l = 'afrac:lfrac:lfrin'
!    character(*),parameter :: fraclist_g = 'gfrac:lfrac'
!    character(*),parameter :: fraclist_r = 'lfrac:lfrin:rfrac'
!
!  we assume ocean and ice are on the same grids, same masks
!  we assume ocn2atm and ice2atm are masked maps
!  we assume lnd2atm is a global map
!  we assume that the ice fraction evolves in time but that
!    the land model fraction does not.  the ocean fraction then
!    is just the complement of the ice fraction over the region
!    of the ocean/ice mask.
!  we assume that component domains are filled with the total
!    potential mask/fraction on that grid, but that the fractions
!    sent at run time are always the relative fraction covered.
!    for example, if an ice cell can be up to 50% covered in
!    ice and 50% land, then the ice domain should have a fraction
!    value of 0.5 at that grid cell.  at run time though, the ice
!    fraction will be between 0.0 and 1.0 meaning that grid cells
!    is covered with between 0.0 and 0.5 by ice.  the "relative" fractions
!    sent at run-time are corrected by the model to be total fractions
!    such that
!  in general, on every grid,
!              fractions_*(afrac) = 1.0
!              fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
!  where fractions_* are a bundle of fractions on a particular grid and
!    *frac (ie afrac) is the fraction of a particular component in the bundle.
!
!  fraclist_g and fraclist_r don't yet interact with atm, lnd, ice, ocn.
!
!  the fractions are computed fundamentally as follows (although the
!    detailed implementation might be slightly different)
!  initialization (frac_init):
!    afrac is set on all grids
!      fractions_a(afrac) = 1.0
!      fractions_o(afrac) = mapa2o(fractions_a(afrac))
!      fractions_i(afrac) = mapa2i(fractions_a(afrac))
!      fractions_l(afrac) = mapa2l(fractions_a(afrac))
!    initially assume ifrac on all grids is zero
!      fractions_*(ifrac) = 0.0
!    fractions/masks provided by surface components
!      fractions_o(ofrac) = dom_o(frac)  ! ocean "mask"
!      fractions_l(lfrin) = dom_l(frac)  ! land model fraction
!    then mapped to the atm model
!      fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
!      fractions_a(lfrin) = mapl2a(fractions_l(lfrin))
!    and a few things are then derived
!      fractions_a(lfrac) = 1.0 - fractions_a(ofrac)
!        this is truncated to zero for very small values (< 0.001)
!        to attempt to preserve non-land gridcells.
!      fractions_l(lfrac) = mapa2l(fractions_a(lfrac))
!      fractions_r(lfrac) = mapl2r(fractions_l(lfrac))
!      fractions_r(lfrin) = mapl2r(fractions_l(lfrin))
!      fractions_g(lfrac) = mapl2g(fractions_l(lfrac))
!
!  run-time (frac_set):
!    update fractions on ice grid
!      fractions_i(ifrac) = i2x_i(Si_ifrac)  ! ice frac from ice model
!      fractions_i(ofrac) = 1.0 - fractions_i(ifrac)
!        note: the relative fractions are corrected to total fractions
!      fractions_o(ifrac) = mapi2o(fractions_i(ifrac))
!      fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
!      fractions_a(ifrac) = mapi2a(fractions_i(ifrac))
!      fractions_a(ofrac) = mapi2a(fractions_i(ofrac))
!
!  fractions used in merging are as follows
!    mrg_x2a uses fractions_a(lfrac,ofrac,ifrac)
!    mrg_x2o needs to use fractions_o(ofrac,ifrac) normalized to one
!      normalization happens in mrg routine
!
!  fraction corrections in mapping are as follows
!    mapo2a uses *fractions_o(ofrac) and /fractions_a(ofrac)
!    mapi2a uses *fractions_i(ifrac) and /fractions_a(ifrac)
!    mapl2a uses *fractions_l(lfrin) and /fractions_a(lfrin)
!    mapl2g weights by fractions_l(lfrac) with normalization, and multiplies by
!      fractions_g(lfrac)
!    mapa2* should use *fractions_a(afrac) and /fractions_*(afrac) but this
!      has been defered since the ratio always close to 1.0
!
!  budgets use the standard afrac, ofrac, ifrac, and lfrac to compute
!
!  NOTE: In trigrid configurations, lfrin MUST be defined as the 
!  conservative o2l mapping of the complement of the ocean mask.
!  In non-trigrid configurations, lfrin is generally associated with
!  the fraction of land grid defined by the surface dataset and might
!  be 1 everywhere for instance.  In many cases, the non-trigrid
!  lfrin is defined to be the conservative o2a mapping of the complement
!  of the ocean mask.  In this case, it is defined the same as the
!  trigrid.  But to support all cases, 
!  for trigrid:
!    mapping from the land grid should use the lfrin field (same in non-trigrid)
!    budget diagnostics should use lfrin (lfrac in non-trigrid)
!    merges in the atm should use lfrac (same in non-trigrid)
!    the runoff should use the lfrin fraction in the runoff merge (lfrac in non-trigrid)
!
!  fraction and domain checks
!    initialization:
!      dom_i = mapo2i(dom_o)  ! lat, lon, mask, area
!      where fractions_a(lfrac) > 0.0, fractions_a(lfrin) is also > 0.0
!         this ensures the land will provide data everywhere the atm needs it
!         and allows the land frac to be subtlely different from the
!         land fraction specified in the atm.
!      dom_a = mapl2a(dom_l)  ! if atm/lnd same grids
!      dom_a = mapo2a(dom_o)  ! if atm/ocn same grids
!      dom_a = mapi2a(dom_i)  ! if atm/ocn same grids
!      0.0-eps < fractions_*(*) < 1.0+eps
!    run time:
!      fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
!      0.0-eps < fractions_*(*) < 1.0+eps
!
!!
!
! !REVISION HISTORY:
!    2007-may-07 - M. Vertenstein - initial port to cpl7.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_frac_mct

  ! !USES:

  use shr_kind_mod     , only: R8 => SHR_KIND_R8
  use shr_sys_mod
  use shr_const_mod

  use mct_mod
  use seq_infodata_mod
  use seq_comm_mct,      only: logunit, loglevel, seq_comm_mpicom, seq_comm_iamroot, CPLID
  use seq_map_mod,       only: seq_map_map
  use seq_map_type_mod,  only: seq_map

  use prep_lnd_mod, only: prep_lnd_get_mapper_Fa2l
  use prep_ocn_mod, only: prep_ocn_get_mapper_Fa2o
  use prep_ocn_mod, only: prep_ocn_get_mapper_SFi2o
  use prep_ice_mod, only: prep_ice_get_mapper_SFo2i
  use prep_rof_mod, only: prep_rof_get_mapper_Fl2r
  use prep_atm_mod, only: prep_atm_get_mapper_Fo2a
  use prep_atm_mod, only: prep_atm_get_mapper_Fi2a
  use prep_atm_mod, only: prep_atm_get_mapper_Fl2a
  use prep_glc_mod, only: prep_glc_get_mapper_Fl2g

  use component_type_mod

  use iMOAB, only: iMOAB_DefineTagStorage
  use seq_comm_mct, only : mbaxid  !           iMOAB id for atm migrated mesh to coupler pes
  use seq_comm_mct, only : mblxid !            iMOAB app id for lnd on cpl pes
#ifdef MOABDEBUG
  use seq_comm_mct, only : mblx2id !           iMOAB id for land mesh instanced from MCT on coupler pes
  use seq_comm_mct, only : mbox2id !           iMOAB id for ocn mesh instanced from MCT on coupler pes
#endif
  use seq_comm_mct, only : mboxid !            iMOAB app id for ocn on cpl pes
  use seq_comm_mct, only : mbofxid !           iMOAB id for mpas ocean migrated mesh to coupler pes, just for xao flux calculations
  use seq_comm_mct, only : mbixid !            iMOAB for sea-ice migrated to coupler
  use seq_comm_mct, only : atm_pg_active !     flag if PG mesh instanced
  use seq_comm_mct, only : mbintxao ! iMOAB id for intx mesh between ocean and atmosphere
  ! for tri grid, sameg_al would be false 

  use seq_comm_mct, only : mbrxid   !          iMOAB id of moab rof migrated to coupler pes 

  use iMOAB, only : iMOAB_DefineTagStorage, iMOAB_SetDoubleTagStorage, &
        iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorageWithGid, iMOAB_WriteMesh, &
        iMOAB_ApplyScalarProjectionWeights, iMOAB_SendElementTag, iMOAB_ReceiveElementTag, &
        iMOAB_FreeSenderBuffers, iMOAB_GetDoubleTagStorage
#ifdef MOABDEBUG
   use seq_comm_mct,     only : num_moab_exports
#endif
  
  use shr_kind_mod, only: CL => SHR_KIND_CL, CX => SHR_KIND_CX, CXX => SHR_KIND_CXX
  use iso_c_binding ! C_NULL_CHAR 
  implicit none
  private
  save

  ! !PUBLIC TYPES:

  ! none

  ! !PUBLIC MEMBER FUNCTIONS:

  public seq_frac_init
  public seq_frac_set

  ! !PUBLIC DATA MEMBERS:

  !EOP

  ! !LOCAL DATA

  integer, private :: seq_frac_debug = 1
  logical, private :: seq_frac_abort = .true.
  logical, private :: seq_frac_dead

  integer :: local_size_mb_ocn ! use it to set/get ofrac from ocn to second ocean instance

  !--- standard ---
  real(r8),parameter :: eps_fracsum = 1.0e-02   ! allowed error in sum of fracs
  real(r8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  real(r8),parameter :: eps_fraclim = 1.0e-03   ! truncation limit in fractions_a(lfrac)
  logical ,parameter :: atm_frac_correct = .false. ! turn on frac correction on atm grid
  !--- standard plus atm fraction consistency ---
  !  real(r8),parameter :: eps_fracsum = 1.0e-12   ! allowed error in sum of fracs
  !  real(r8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  !  real(r8),parameter :: eps_fraclim = 1.0e-03   ! truncation limit in fractions_a(lfrac)
  !  logical ,parameter :: atm_frac_correct = .true. ! turn on frac correction on atm grid
  !--- unconstrained and area conserving? ---
  !  real(r8),parameter :: eps_fracsum = 1.0e-12   ! allowed error in sum of fracs
  !  real(r8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  !  real(r8),parameter :: eps_fraclim = 1.0e-20   ! truncation limit in fractions_a(lfrac)
  !  logical ,parameter :: atm_frac_correct = .true. ! turn on frac correction on atm grid

  type(seq_map)  , pointer :: mapper_o2a
  type(seq_map)  , pointer :: mapper_i2a
  type(seq_map)  , pointer :: mapper_l2a
  type(seq_map)  , pointer :: mapper_o2i
  type(seq_map)  , pointer :: mapper_a2o
  type(seq_map)  , pointer :: mapper_i2o
  type(seq_map)  , pointer :: mapper_a2l
  type(seq_map)  , pointer :: mapper_l2r
  type(seq_map)  , pointer :: mapper_l2g

  private seq_frac_check

  !===============================================================================
contains
  !===============================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_frac_init
  !
  ! !DESCRIPTION:
  !    Initialize fraction attribute vectors and necessary ocn/ice domain
  !    fraction input if appropriate
  !
  ! !REVISION HISTORY:
  !    2007-may-07 - M. Vertenstein - initial cpl7 version.
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_frac_init( infodata,         &
       atm, ice, lnd, ocn, glc, rof, wav, iac,&
       fractions_a, fractions_i, fractions_l, &
       fractions_o, fractions_g, fractions_r, &
       fractions_w, fractions_z)

    ! !INPUT/OUTPUT PARAMETERS:
    type(seq_infodata_type) , intent(in)    :: infodata
    type(component_type)    , intent(in)    :: atm
    type(component_type)    , intent(in)    :: ice
    type(component_type)    , intent(in)    :: lnd
    type(component_type)    , intent(in)    :: ocn
    type(component_type)    , intent(in)    :: glc
    type(component_type)    , intent(in)    :: rof
    type(component_type)    , intent(in)    :: wav
    type(component_type)    , intent(in)    :: iac
    type(mct_aVect)         , intent(inout) :: fractions_a   ! Fractions on atm grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_i   ! Fractions on ice grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_l   ! Fractions on lnd grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_o   ! Fractions on ocn grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_g   ! Fractions on glc grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_r   ! Fractions on rof grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_w   ! Fractions on wav grid/decomp
    type(mct_aVect)         , intent(inout) :: fractions_z   ! Fractions on iac grid/decomp
    !EOP

    !----- local -----
    type(mct_ggrid), pointer    :: dom_a
    type(mct_ggrid), pointer    :: dom_i
    type(mct_ggrid), pointer    :: dom_l
    type(mct_ggrid), pointer    :: dom_o
    type(mct_ggrid), pointer    :: dom_g
    type(mct_ggrid), pointer    :: dom_r
    type(mct_ggrid), pointer    :: dom_w
    type(mct_ggrid), pointer    :: dom_z

    logical :: atm_present   ! .true. => atm is present
    logical :: ice_present   ! .true. => ice is present
    logical :: ocn_present   ! .true. => ocean is present
    logical :: lnd_present   ! .true. => land is present
    logical :: glc_present   ! .true. => glc is present
    logical :: rof_present   ! .true. => rof is present
    logical :: wav_present   ! .true. => wav is present
    logical :: iac_present   ! .true. => iac is present
    logical :: dead_comps    ! .true. => dead models present

    integer :: n            ! indices
    integer :: ka, ki, kl, ko ! indices
    integer :: kf, kk, kr, kg ! indices
    integer :: lsize          ! local size of ice av
    integer :: debug_old      ! old debug value

    character(*),parameter :: fraclist_a = 'afrac:ifrac:ofrac:lfrac:lfrin'
    character(*),parameter :: fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
    character(*),parameter :: fraclist_i = 'afrac:ifrac:ofrac'
    character(*),parameter :: fraclist_l = 'afrac:lfrac:lfrin'
    character(*),parameter :: fraclist_g = 'gfrac:lfrac'
    character(*),parameter :: fraclist_r = 'lfrac:lfrin:rfrac'
    character(*),parameter :: fraclist_w = 'wfrac'
    character(*),parameter :: fraclist_z = 'afrac:lfrac'

    ! moab
    integer                  :: tagtype, numco,  tagindex, ent_type, ierr, arrSize
    character(CXX)           :: tagname, tagNameExt
    character*32             ::  wgtIdef
    real(r8),    allocatable :: tagValues(:) ! used for setting some default tags
    integer ,    allocatable :: GlobalIds(:) ! used for setting values associated with ids
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
    integer kgg  ! index in global number attribute, used for global id in MOAB 
    integer idintx ! used for context for intx atm - ocn
    integer id_join ! used for example for atm%cplcompid
    integer :: mpicom ! we are on coupler PES here
    character(30)            :: outfile, wopts


    !----- formats -----
    character(*),parameter :: subName = '(seq_frac_init) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         atm_present=atm_present,       &
         lnd_present=lnd_present,       &
         rof_present=rof_present,       &
         ice_present=ice_present,       &
         ocn_present=ocn_present,       &
         glc_present=glc_present,       &
         wav_present=wav_present,       &
         iac_present=iac_present,       &
         dead_comps=dead_comps)

    dom_a => component_get_dom_cx(atm)
    dom_l => component_get_dom_cx(lnd)
    dom_i => component_get_dom_cx(ice)
    dom_o => component_get_dom_cx(ocn)
    dom_r => component_get_dom_cx(rof)
    dom_g => component_get_dom_cx(glc)
    dom_w => component_get_dom_cx(wav)
    dom_z => component_get_dom_cx(iac)

    debug_old = seq_frac_debug
    seq_frac_debug = 2

    ! Initialize fractions on atm grid/decomp (initialize ice fraction to zero)

    if (atm_present) then
       lSize = mct_aVect_lSize(dom_a%data)
       call mct_aVect_init(fractions_a,rList=fraclist_a,lsize=lsize)
       call mct_aVect_zero(fractions_a)

       ka = mct_aVect_indexRa(fractions_a,"afrac",perrWith=subName)
       fractions_a%rAttr(ka,:) = 1.0_r8

       ! Initialize fractions on atm coupler mesh; on migrated atm to coupler
       if (mbaxid .ge. 0  ) then ! // 
         tagname = trim(fraclist_a)//C_NULL_CHAR ! 'afrac:ifrac:ofrac:lfrac:lfrin'
         tagtype = 1  ! dense, double
         numco = 1 !   
         ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags on atm phys mesh on cpl '
            call shr_sys_abort(subname//' ERROR in defining tags on atm phys mesh on cpl')
         endif
         ! find out the number of local elements in moab mesh
         ierr  = iMOAB_GetMeshInfo ( mbaxid, nvert, nvise, nbl, nsurf, nvisBC );
         
         ! we should set to 1 the 'afrac' tag
         ! we are on cells now !
         arrSize = nvise(1) * 5 ! there are 5 tags that need to be zeroed out 
         allocate(tagValues(arrSize) )
         ent_type = 1 ! cell type
         tagValues = 0.0_r8
         ierr = iMOAB_SetDoubleTagStorage ( mbaxid, tagname, arrSize , ent_type, tagValues)
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in zeroing out fracs  '
            call shr_sys_abort(subname//' ERROR in zeroing out fracs on phys atm')
         endif
         deallocate(tagValues)


         allocate(tagValues(nvise(1)))
         tagname = 'afrac'//C_NULL_CHAR
         tagValues = 1.0_r8
         ierr = iMOAB_SetDoubleTagStorage ( mbaxid, tagname, nvise(1) , ent_type, tagValues)


         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in setting afrac tag on phys atm  '
            call shr_sys_abort(subname//' ERROR in setting afrac tag on phys atm')
         endif
         deallocate(tagValues)
       endif
    endif

    ! Initialize fractions on glc grid decomp, just an initial "guess", updated later

    if (glc_present) then
       lSize = mct_aVect_lSize(dom_g%data)
       call mct_aVect_init(fractions_g,rList=fraclist_g,lsize=lsize)
       call mct_aVect_zero(fractions_g)

       kg = mct_aVect_indexRA(fractions_g,"gfrac",perrWith=subName)
       kf = mct_aVect_indexRA(dom_g%data ,"frac" ,perrWith=subName)
       fractions_g%rAttr(kg,:) = dom_g%data%rAttr(kf,:)
    end if

    ! Initialize fractions on land grid decomp, just an initial "guess", updated later

    if (lnd_present) then
       lSize = mct_aVect_lSize(dom_l%data)
       call mct_aVect_init(fractions_l,rList=fraclist_l,lsize=lsize)
       call mct_aVect_zero(fractions_l)
       if (mblxid .ge. 0  ) then ! // 
         tagname = trim(fraclist_l)//C_NULL_CHAR ! 'afrac:lfrac:lfrin'
         tagtype = 1  ! dense, double
         numco = 1 !   
         ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags on lnd phys mesh on cpl '
            call shr_sys_abort(subname//' ERROR in defining tags on lnd phys mesh on cpl')
         endif
         ierr  = iMOAB_GetMeshInfo ( mblxid, nvert, nvise, nbl, nsurf, nvisBC );
         ! this should have some cells 
         ent_type = 1 ! cells from now on
         arrSize = 3 * nvise(1) 

         allocate(tagValues(arrSize) )
         tagValues = 0 
         ierr = iMOAB_SetDoubleTagStorage ( mblxid, tagname, arrSize , ent_type, tagValues)

         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in setting fractions tags on lnd   '
            call shr_sys_abort(subname//' ERROR in setting fractions tags on lnd')
         endif
         deallocate(tagValues)
       endif
#ifdef MOABDEBUG
       ! mblx2id is the id for moab app exposing land cpl 
       call expose_mct_grid_moab(lnd, mblx2id)
#endif

       kk = mct_aVect_indexRA(fractions_l,"lfrin",perrWith=subName)
       kf = mct_aVect_indexRA(dom_l%data ,"frac" ,perrWith=subName)
       fractions_l%rAttr(kk,:) = dom_l%data%rAttr(kf,:)
       if (mblxid .ge. 0  ) then ! // 
         tagname = 'lfrin'//C_NULL_CHAR ! 'lfrin'
         allocate(tagValues(lSize) )
         tagValues = dom_l%data%rAttr(kf,:)
         kgg = mct_aVect_indexIA(dom_l%data ,"GlobGridNum" ,perrWith=subName)
         allocate(GlobalIds(lSize))
         GlobalIds = dom_l%data%iAttr(kgg,:)

         ! ent_type should be 3, FV
         ierr = iMOAB_SetDoubleTagStorageWithGid ( mblxid, tagname, lSize , ent_type, tagValues, GlobalIds )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in setting lfrin on lnd   '
            call shr_sys_abort(subname//' ERROR in setting lfrin on lnd')
         endif
         deallocate(GlobalIds)
         deallocate(tagValues)

       endif

       if (atm_present) then
          mapper_l2a => prep_atm_get_mapper_Fl2a()
          mapper_a2l => prep_lnd_get_mapper_Fa2l()
          call seq_map_map(mapper_l2a, fractions_l, fractions_a, fldlist='lfrin', norm=.false.)
          call seq_map_map(mapper_a2l, fractions_a, fractions_l, fldlist='afrac', norm=.false.)
       endif

    end if  ! end of (if lnd_present)

    ! Initialize fractions on ice grid/decomp (initialize ice fraction to zero)

    if (rof_present) then
       lSize = mct_aVect_lSize(dom_r%data)
       call mct_aVect_init(fractions_r,rList=fraclist_r,lsize=lsize)
       call mct_aVect_zero(fractions_r)

       kr = mct_aVect_indexRa(fractions_r,"rfrac",perrWith=subName)
       kf = mct_aVect_indexRA(dom_r%data ,"frac" ,perrWith=subName)
       fractions_r%rAttr(kr,:) = dom_r%data%rAttr(kf,:)
       ! fractions needs to be set on river grid (that has a mask ?)
       ! copy from land code; we have on coupler the river instance that comes from river-ocean map
       if (mbrxid .ge. 0  ) then ! // 
         ! set all to zero, and then modify rfrac
         tagname = trim(fraclist_r)//C_NULL_CHAR ! 'lfrac:lfrin:rfrac'
         tagtype = 1  ! dense, double
         numco = 1 !   
         ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags on rof phys mesh on cpl '
            call shr_sys_abort(subname//' ERROR in defining tags on rof phys mesh on cpl')
         endif
         ierr  = iMOAB_GetMeshInfo ( mbrxid, nvert, nvise, nbl, nsurf, nvisBC );
         arrSize = 3 * nvise(1) ! there are 3 tags 
         allocate(tagValues(arrSize) )
         ent_type = 1 ! vertex type, rof is now FV
         tagValues = 0. 
         ierr = iMOAB_SetDoubleTagStorage ( mbrxid, tagname, arrSize , ent_type, tagValues)
         deallocate(tagValues)

         tagname = 'rfrac'//C_NULL_CHAR ! 'rfrac'
         allocate(tagValues(lSize) )
         tagValues = dom_r%data%rAttr(kf,:)
         kgg = mct_aVect_indexIA(dom_r%data ,"GlobGridNum" ,perrWith=subName)
         allocate(GlobalIds(lSize))
         GlobalIds = dom_r%data%iAttr(kgg,:)
         ! again, we are setting on the river instance that is also used for ocean coupling
         ierr = iMOAB_SetDoubleTagStorageWithGid ( mbrxid, tagname, lSize , ent_type, tagValues, GlobalIds )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in setting rfrac on rof   '
            call shr_sys_abort(subname//' ERROR in setting rfrac on rof ')
         endif
         deallocate(GlobalIds)
         deallocate(tagValues)

       endif
    end if

    ! Initialize fractions on wav grid decomp, just an initial "guess", updated later

    if (wav_present) then
       lSize = mct_aVect_lSize(dom_w%data)
       call mct_aVect_init(fractions_w,rList=fraclist_w,lsize=lsize)
       call mct_aVect_zero(fractions_w)
       fractions_w%rAttr(:,:) = 1.0_r8
    end if

    ! Initialize fractions on iac grid decomp, just an initial "guess", updated later

    if (iac_present) then
       lSize = mct_aVect_lSize(dom_z%data)
       call mct_aVect_init(fractions_z,rList=fraclist_z,lsize=lsize)
       call mct_aVect_zero(fractions_z)
       fractions_z%rAttr(:,:) = 1.0_r8
    end if

    ! Initialize fractions on ice grid/decomp (initialize ice fraction to zero)

    if (ice_present) then
       lSize = mct_aVect_lSize(dom_i%data)
       call mct_aVect_init(fractions_i,rList=fraclist_i,lsize=lsize)
       call mct_aVect_zero(fractions_i)

       ko = mct_aVect_indexRa(fractions_i,"ofrac",perrWith=subName)
       kf = mct_aVect_indexRA(dom_i%data ,"frac" ,perrWith=subName)
       fractions_i%rAttr(ko,:) = dom_i%data%rAttr(kf,:)

       if (mbixid .ge. 0  ) then ! // 
         tagname = trim(fraclist_i)//C_NULL_CHAR ! 'afrac:ifrac:ofrac'
         tagtype = 1  ! dense, double
         numco = 1 !   
         ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags on ice phys mesh on cpl '
            call shr_sys_abort(subname//' ERROR in defining tags on ice phys mesh on cpl')
         endif
         ! start copy from rof
         ierr  = iMOAB_GetMeshInfo ( mbixid, nvert, nvise, nbl, nsurf, nvisBC );
         arrSize = 3 * nvise(1) ! there are 3 tags 'afrac:ifrac:ofrac'
         allocate(tagValues(arrSize) )
         ent_type = 1  ! cell type, ice is FV
         tagValues = 0. 
         ierr = iMOAB_SetDoubleTagStorage ( mbixid, tagname, arrSize , ent_type, tagValues)
         deallocate(tagValues)
         tagname = 'ofrac'//C_NULL_CHAR ! 'ofrac'
         allocate(tagValues(lSize) )
         tagValues = dom_i%data%rAttr(kf,:)
         kgg = mct_aVect_indexIA(dom_i%data ,"GlobGridNum" ,perrWith=subName)
         allocate(GlobalIds(lSize))
         GlobalIds = dom_i%data%iAttr(kgg,:)

         ierr = iMOAB_SetDoubleTagStorageWithGid ( mbixid, tagname, lSize , ent_type, tagValues, GlobalIds )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in setting ofrac on ice   '
            call shr_sys_abort(subname//' ERROR in setting ofrac on ice ')
         endif
         deallocate(GlobalIds)
         deallocate(tagValues)
       endif

       if (atm_present) then
          mapper_i2a => prep_atm_get_mapper_Fi2a()
          call seq_map_map(mapper_i2a,fractions_i,fractions_a,fldlist='ofrac',norm=.false.)
       endif

    end if ! end of ice_present

    ! Initialize fractions on ocean grid/decomp (initialize ice fraction to zero)
    ! These are initialized the same as for ice

    if (ocn_present) then
       lSize = mct_aVect_lSize(dom_o%data)
       call mct_aVect_init(fractions_o,rList=fraclist_o,lsize=lsize)
       call mct_aVect_zero(fractions_o)
#ifdef MOABDEBUG
       ! initialize ocn imoab app on mct grid 
       call expose_mct_grid_moab(ocn, mbox2id) ! will use then to set the data on it , for debugging
#endif
       if (mboxid .ge. 0  ) then ! // 
         tagname = trim(fraclist_o)//C_NULL_CHAR ! 'afrac:ifrac:ofrac:ifrad:ofrad'
         tagtype = 1  ! dense, double
         numco = 1 !   
         ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags on ocn phys mesh on cpl '
            call shr_sys_abort(subname//' ERROR in defining tags on ocn phys mesh on cpl')
         endif
         ierr = iMOAB_DefineTagStorage(mbofxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags on ocn phys mesh on cpl '
            call shr_sys_abort(subname//' ERROR in defining tags on ocn phys mesh on cpl')
         endif
          ierr  = iMOAB_GetMeshInfo ( mboxid, nvert, nvise, nbl, nsurf, nvisBC );
         arrSize = 5 * nvise(1) ! there are 5 tags 'afrac:ifrac:ofrac:ifrad:ofrad'
         local_size_mb_ocn = nvise(1)
         allocate(tagValues(arrSize) )
         ent_type = 1  ! cell type, ocn is FV
         tagValues = 0. 
         ierr = iMOAB_SetDoubleTagStorage ( mboxid, tagname, arrSize , ent_type, tagValues)
         ierr = iMOAB_SetDoubleTagStorage ( mbofxid, tagname, arrSize , ent_type, tagValues)
         deallocate(tagValues)

       endif
       if (ice_present) then
          mapper_i2o => prep_ocn_get_mapper_SFi2o()
          call seq_map_map(mapper_i2o,fractions_i,fractions_o,fldlist='ofrac',norm=.false.)

          tagname = 'ofrac'//C_NULL_CHAR
          ent_type = 1! cells
          allocate(tagValues(local_size_mb_ocn) )
          ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, local_size_mb_ocn , ent_type, tagValues)
          ierr = iMOAB_SetDoubleTagStorage ( mbofxid, tagname, local_size_mb_ocn , ent_type, tagValues)
          deallocate(tagValues)

       else
       ! still need to TODO moab case when no sea ice
          ko = mct_aVect_indexRa(fractions_o,"ofrac",perrWith=subName)
          kf = mct_aVect_indexRA(dom_o%data ,"frac" ,perrWith=subName)
          fractions_o%rAttr(ko,:) = dom_o%data%rAttr(kf,:)
          ! so the frac var from dom is copied into ofrac field
          !  frac is on ocean mesh, ofrac is on ocean mesh too; 
          tagname = 'frac'//C_NULL_CHAR
          ent_type = 1! cells
          allocate(tagValues(local_size_mb_ocn) )
          ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, local_size_mb_ocn , ent_type, tagValues)
          tagname = 'ofrac'//C_NULL_CHAR
          ierr = iMOAB_SetDoubleTagStorage ( mboxid, tagname, local_size_mb_ocn , ent_type, tagValues)
          ierr = iMOAB_SetDoubleTagStorage ( mbofxid, tagname, local_size_mb_ocn , ent_type, tagValues)
          deallocate(tagValues)
          mapper_o2a => prep_atm_get_mapper_Fo2a()
          call seq_map_map(mapper_o2a, fractions_o, fractions_a, fldlist='ofrac',norm=.false.)
       endif

       if (atm_present) then
          mapper_a2o => prep_ocn_get_mapper_Fa2o()
          call seq_map_map(mapper_a2o, fractions_a, fractions_o, fldlist='afrac',norm=.false.)
       endif


       if (ice_present) then
          ! --- this should be an atm2ice call above, but atm2ice doesn't work
          mapper_o2i => prep_ice_get_mapper_SFo2i()
          call seq_map_map(mapper_o2i,fractions_o,fractions_i,fldlist='afrac',norm=.false.)
       endif
    end if  ! end of if ocn present

    ! --- Set ofrac and lfrac on atm grid.  These should actually be mapo2a of
    !     ofrac and lfrac but we can't map lfrac from o2a due to masked mapping
    !     weights.  So we have to settle for a residual calculation that is
    !     truncated to zero to try to preserve "all ocean" cells.

    if (atm_present) then
       ka = mct_aVect_indexRa(fractions_a,"afrac",perrWith=subName)
       ki = mct_aVect_indexRa(fractions_a,"ifrac",perrWith=subName)
       kl = mct_aVect_indexRa(fractions_a,"lfrac",perrWith=subName)
       ko = mct_aVect_indexRa(fractions_a,"ofrac",perrWith=subName)
       kk = mct_aVect_indexRa(fractions_a,"lfrin",perrWith=subName)
       lSize = mct_aVect_lSize(fractions_a)

       if (ice_present .or. ocn_present) then
          do n = 1,lsize
             fractions_a%rAttr(kl,n) = 1.0_r8 - fractions_a%rAttr(ko,n)
             if (abs(fractions_a%rAttr(kl,n)) < eps_fraclim) then
                fractions_a%rAttr(kl,n) = 0.0_r8
                if (atm_frac_correct) fractions_a%rAttr(ko,n) = 1.0_r8
             endif
          enddo
       ! TODO: replace this with math
        if (mbaxid .ge. 0  ) then ! //
          tagname = 'lfrac'//C_NULL_CHAR ! 'lfrac
          allocate(tagValues(lSize) )
          tagValues = fractions_a%rAttr(kl,:)
          kgg = mct_aVect_indexIA(dom_a%data ,"GlobGridNum" ,perrWith=subName)
          allocate(GlobalIds(lSize))
          GlobalIds = dom_a%data%iAttr(kgg,:)
          ! set on atmosphere instance
          ierr = iMOAB_SetDoubleTagStorageWithGid ( mbaxid, tagname, lSize , ent_type, tagValues, GlobalIds )
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in setting lfrac on atm   '
             call shr_sys_abort(subname//' ERROR in setting lfrac on atm ')
          endif
          deallocate(GlobalIds)
          deallocate(tagValues)
        endif
       else if (lnd_present) then
       ! TODO:  MOAB case.
          do n = 1,lsize
             fractions_a%rAttr(kl,n) = fractions_a%rAttr(kk,n)
             fractions_a%rAttr(ko,n) = 1.0_r8 - fractions_a%rAttr(kl,n)
             if (abs(fractions_a%rAttr(ko,n)) < eps_fraclim) then
                fractions_a%rAttr(ko,n) = 0.0_r8
                if (atm_frac_correct) fractions_a%rAttr(kl,n) = 1.0_r8
             endif
          enddo
       endif
    endif

    ! --- finally, set fractions_l(lfrac) from fractions_a(lfrac)
    ! --- and fractions_r(lfrac:lfrin) from fractions_l(lfrac:lfrin)
    ! --- and fractions_g(lfrac) from fractions_l(lfrac)

    if (lnd_present) then
       if (atm_present) then
          mapper_a2l => prep_lnd_get_mapper_Fa2l()
          call seq_map_map(mapper_a2l, fractions_a, fractions_l, fldlist='lfrac', norm=.false.)
       else
          ! If the atmosphere is absent, then simply set fractions_l(lfrac) = fractions_l(lfrin)
          kk = mct_aVect_indexRA(fractions_l,"lfrin",perrWith=subName)
          kl = mct_aVect_indexRA(fractions_l,"lfrac",perrWith=subName)
          fractions_l%rAttr(kl,:) = fractions_l%rAttr(kk,:)
       end if
    end if
    if (lnd_present .and. rof_present) then
       mapper_l2r => prep_rof_get_mapper_Fl2r()
       call seq_map_map(mapper_l2r, fractions_l, fractions_r, fldlist='lfrac:lfrin', norm=.false.)
    endif
    if (lnd_present .and. glc_present) then
       mapper_l2g => prep_glc_get_mapper_Fl2g()
       call seq_map_map(mapper_l2g, fractions_l, fractions_g, fldlist='lfrac', norm=.false.)
    end if

    if (lnd_present) call seq_frac_check(fractions_l,'lnd init')
    if (glc_present) call seq_frac_check(fractions_g,'glc init')
    if (rof_present) call seq_frac_check(fractions_r,'rof init')
    if (wav_present) call seq_frac_check(fractions_w,'wav init')
    if (iac_present) call seq_frac_check(fractions_z,'iac init')
    if (ice_present) call seq_frac_check(fractions_i,'ice init')
    if (ocn_present) call seq_frac_check(fractions_o,'ocn init')
    if (atm_present .and. (lnd_present.or.ice_present.or.ocn_present)) &
         call seq_frac_check(fractions_a,'atm init')
    seq_frac_debug = debug_old
#ifdef MOABDEBUG
    wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR
    if (mbaxid .ge. 0  ) then
        outfile = 'atmCplFr.h5m'//C_NULL_CHAR
        ierr = iMOAB_WriteMesh(mbaxid, trim(outfile), trim(wopts))
        if (ierr .ne. 0) then
          write(logunit,*) subname,' error in writing mesh '
          call shr_sys_abort(subname//' ERROR in writing mesh ')
        endif
    endif
    if (mblxid .ge. 0  ) then
        outfile = 'lndCplFr.h5m'//C_NULL_CHAR
        ierr = iMOAB_WriteMesh(mblxid, trim(outfile), trim(wopts))
        if (ierr .ne. 0) then
          write(logunit,*) subname,' error in writing mesh '
          call shr_sys_abort(subname//' ERROR in writing mesh ')
        endif
    endif
    if (mbrxid .ge. 0  ) then
        outfile = 'rofCplFr.h5m'//C_NULL_CHAR
        ierr = iMOAB_WriteMesh(mbrxid, trim(outfile), trim(wopts))
        if (ierr .ne. 0) then
          write(logunit,*) subname,' error in writing mesh '
          call shr_sys_abort(subname//' ERROR in writing mesh ')
        endif
    endif
    if (mbixid .ge. 0  ) then
        outfile = 'iceCplFr.h5m'//C_NULL_CHAR
        ierr = iMOAB_WriteMesh(mbixid, trim(outfile), trim(wopts))
        if (ierr .ne. 0) then
          write(logunit,*) subname,' error in writing mesh '
          call shr_sys_abort(subname//' ERROR in writing mesh ')
        endif
    endif
    if (mboxid .ge. 0  ) then
        outfile = 'ocnCplFr.h5m'//C_NULL_CHAR
        ierr = iMOAB_WriteMesh(mboxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing mesh '
            call shr_sys_abort(subname//' ERROR in writing mesh ')
         endif
    endif
#endif

  end subroutine seq_frac_init

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_frac_set
  !
  ! !DESCRIPTION:
  !    Update surface fractions
  !
  ! !REVISION HISTORY:
  !    2007-may-07 - M. Vertenstein - initial cpl7 version.
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_frac_set(infodata, ice, ocn, fractions_a, fractions_i, fractions_o)

    use seq_comm_mct, only : mboxid !            iMOAB app id for ocn on cpl pes

    ! !INPUT/OUTPUT PARAMETERS:
    type(seq_infodata_type) , intent(in)    :: infodata
    type(component_type)    , intent(in)    :: ice
    type(component_type)    , intent(in)    :: ocn
    type(mct_aVect)         , intent(inout) :: fractions_a   ! Fractions on atm
    type(mct_aVect)         , intent(inout) :: fractions_i   ! Fractions on ice
    type(mct_aVect)         , intent(inout) :: fractions_o   ! Fractions on ocn
    !EOP

    !----- local -----
    type(mct_aVect), pointer :: i2x_i
    type(mct_ggrid), pointer :: dom_i
    type(mct_ggrid), pointer :: dom_o ! introduced just to update ocean moab fractions
    logical                  :: atm_present   ! true => atm is present
    logical                  :: ice_present   ! true => ice is present
    logical                  :: ocn_present   ! true => ocn is present
    integer                  :: n
    integer                  :: ki, kl, ko, kf
    real(r8),allocatable :: fcorr(:)

    logical, save :: first_time = .true.
! moab
    integer                  ::   ierr, kgg
    integer , save           ::   lSize, ent_type
    character(CXX)           :: tagname
    real(r8),    allocatable, save :: tagValues(:) ! used for setting some tags
    real(r8),    allocatable, save :: tagValuesOfrac(:) ! used for setting some tags
    integer ,    allocatable, save :: GlobalIds(:) ! used for setting values associated with ids
#ifdef MOABDEBUG
    character(len=100) :: outfile, wopts, lnum
#endif

    !----- formats -----
    character(*),parameter :: subName = '(seq_frac_set) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !----------------------------------------------------------------------
    ! Update fractions
    ! - Update ice fraction on ice grid first, normalize to total fraction
    !   available for cover
    ! - Update ocn fraction on ice grid as residual
    ! - Map ice/ocn fractions from ice grid to ocean and atm grids
    !----------------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         atm_present=atm_present,       &
         ice_present=ice_present,       &
         ocn_present=ocn_present)

    dom_i => component_get_dom_cx(ice)
    i2x_i => component_get_c2x_cx(ice)

    dom_o => component_get_dom_cx(ocn) ! 
    if (ice_present) then
       call mct_aVect_copy(i2x_i, fractions_i, "Si_ifrac", "ifrac")

       ki = mct_aVect_indexRA(fractions_i,"ifrac")
       ko = mct_aVect_indexRA(fractions_i,"ofrac")
       kf = mct_aVect_indexRA(dom_i%data ,"frac" ,perrWith=subName)
       fractions_i%rAttr(ki,:) = fractions_i%rAttr(ki,:) * dom_i%data%rAttr(kf,:)
       fractions_i%rAttr(ko,:) = dom_i%data%rAttr(kf,:) - fractions_i%rAttr(ki,:)

       call seq_frac_check(fractions_i,'ice set')

       ! update ice fractions on moab instance
       if (first_time) then  ! allocate some local arrays
          lSize = mct_aVect_lSize(dom_i%data)
          allocate(tagValues(lSize) )
          allocate(GlobalIds(lSize) )
          kgg = mct_aVect_indexIA(dom_o%data ,"GlobGridNum" ,perrWith=subName)
          GlobalIds = dom_i%data%iAttr(kgg,:)
          allocate (tagValuesOfrac(local_size_mb_ocn))
          ent_type = 1 ! cells for mpas sea ice
          first_time = .false.
       endif

          ! something like this:
       if (mbixid > 0 ) then

          tagname = 'ifrac'//C_NULL_CHAR
          ! fraclist_i = 'afrac:ifrac:ofrac'
          !
          tagValues = fractions_i%rAttr(ki,:)
          ierr = iMOAB_SetDoubleTagStorageWithGid ( mbixid, tagname, lSize , ent_type, tagValues, GlobalIds )
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in setting ifrac on ice moab instance  '
             call shr_sys_abort(subname//' ERROR in setting ifrac on ice moab instance ')
          endif

          tagname = 'ofrac'//C_NULL_CHAR
          tagValues = fractions_i%rAttr(ko,:)
          ierr = iMOAB_SetDoubleTagStorageWithGid ( mbixid, tagname, lSize , ent_type, tagValues, GlobalIds )
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in setting ofrac on ice moab instance  '
             call shr_sys_abort(subname//' ERROR in setting ofrac on ice moab instance ')
          endif
#ifdef MOABDEBUG
         write(lnum,"(I0.2)")num_moab_exports
         outfile = 'iceCplFr_'//trim(lnum)//'.h5m'//C_NULL_CHAR
         wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
         ierr = iMOAB_WriteMesh(mbixid, outfile, wopts)
#endif

       endif

       if (ocn_present) then
          mapper_i2o => prep_ocn_get_mapper_SFi2o()
          call seq_map_map(mapper_i2o, fractions_i, fractions_o, &
               fldlist='ofrac:ifrac',norm=.false.)
          call seq_frac_check(fractions_o, 'ocn set')
          ! set the ofrac artificially on mbofxid instance, because it is needed for prep_aoxflux
          tagname = 'ofrac'//C_NULL_CHAR
          ent_type = 1! cells
          ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, local_size_mb_ocn , ent_type, tagValuesOfrac)
          ierr = iMOAB_SetDoubleTagStorage ( mbofxid, tagname, local_size_mb_ocn , ent_type, tagValuesOfrac)
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in setting ofrac on mbofxid moab instance  '
             call shr_sys_abort(subname//' ERROR in setting ofrac on mbofxid moab instance ')
          endif
       endif


       if (atm_present) then
          mapper_i2a => prep_atm_get_mapper_Fi2a()
          call seq_map_map(mapper_i2a, fractions_i, fractions_a, &
               fldlist='ofrac:ifrac', norm=.false.)
#ifdef MOABDEBUG
            write(lnum,"(I0.2)")num_moab_exports
            outfile = 'atmCplFr_'//trim(lnum)//'.h5m'//C_NULL_CHAR
            wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
            ierr = iMOAB_WriteMesh(mbaxid, outfile, wopts)
#endif
          !tcx---  fraction correction, this forces the fractions_a to sum to 1.0_r8.
          !   ---  but it introduces a conservation error in mapping
          if (atm_frac_correct) then
             ki = mct_aVect_indexRA(fractions_a,"ifrac")
             ko = mct_aVect_indexRA(fractions_a,"ofrac")
             kl = mct_aVect_indexRA(fractions_a,"lfrac")
             lSize = mct_aVect_lSize(fractions_a)
             allocate(fcorr(lsize))
             do n = 1,lsize
                if ((fractions_a%rAttr(ki,n)+fractions_a%rAttr(ko,n)) > 0.0_r8) then
                   fcorr(n) = ((1.0_r8-fractions_a%rAttr(kl,n))/ &
                        (fractions_a%rAttr(ki,n)+fractions_a%rAttr(ko,n)))
                else
                   fcorr(n) = 0.0_r8
                endif
             enddo
             fractions_a%rAttr(ki,:) = fractions_a%rAttr(ki,:) * fcorr(:)
             fractions_a%rAttr(ko,:) = fractions_a%rAttr(ko,:) * fcorr(:)
             deallocate(fcorr)
          endif

          call seq_frac_check(fractions_a,'atm set')
       endif
    end if

  end subroutine seq_frac_set
  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_frac_check
  !
  ! !DESCRIPTION:
  !    Check surface fractions
  !
  ! !REVISION HISTORY:
  !    2008-jun-11 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_frac_check(fractions,string)

    ! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect) , intent(in) :: fractions   ! Fractions datatype
    character(len=*), intent(in), optional :: string      ! character string

    !EOP

    !----- local -----
    integer  :: n, lsize
    integer  :: ncnt
    integer  :: mpicom
    logical  :: iamroot
    real(r8) :: sum,diff,maxerr
    real(r8) :: aminval,amaxval   ! used for lnd
    real(r8) :: lminval,lmaxval   ! used for lnd
    real(r8) :: ominval,omaxval   ! used for ocn
    real(r8) :: iminval,imaxval   ! used for ice
    real(r8) :: gminval,gmaxval   ! used for glc
    real(r8) :: rminval,rmaxval   ! used for rof
    real(r8) :: wminval,wmaxval   ! used for wav
    real(r8) :: zminval,zmaxval   ! used for iac
    real(r8) :: kminval,kmaxval   ! used for lnd, lfrin
    real(r8) :: sminval,smaxval   ! used for sum
    real(r8) :: tmpmin, tmpmax    ! global tmps
    integer  :: tmpsum            ! global tmp
    integer  :: ka,kl,ki,ko,kg,kk,kr,kw,kz
    character(len=128) :: lstring
    logical :: error

    !----- formats -----
    character(*),parameter :: subName = '(seq_frac_check) '
    character(*),parameter :: F01 = "('(seq_frac_check) ',2a,i8,g26.18)"
    character(*),parameter :: F02 = "('(seq_frac_check) ',2a,2g26.18)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    mpicom = seq_comm_mpicom(CPLID)
    iamroot = seq_comm_iamroot(CPLID)

    if (present(string)) then
       lstring='['//trim(string)//']'
    else
       lstring=''
    endif

    ka = -1
    kl = -1
    ki = -1
    ko = -1
    kk = -1
    kg = -1
    kr = -1
    kw = -1
    kz = -1
    aminval =  999.0_r8
    amaxval = -999.0_r8
    lminval =  999.0_r8
    lmaxval = -999.0_r8
    iminval =  999.0_r8
    imaxval = -999.0_r8
    ominval =  999.0_r8
    omaxval = -999.0_r8
    gminval =  999.0_r8
    gmaxval = -999.0_r8
    kminval =  999.0_r8
    kmaxval = -999.0_r8
    sminval =  999.0_r8
    smaxval = -999.0_r8
    rminval =  999.0_r8
    rmaxval = -999.0_r8
    wminval =  999.0_r8
    wmaxval = -999.0_r8
    zminval =  999.0_r8
    zmaxval = -999.0_r8

    lsize = mct_avect_lsize(fractions)
    ka = mct_aVect_indexRA(fractions,"afrac",perrWith='quiet')
    kl = mct_aVect_indexRA(fractions,"lfrac",perrWith='quiet')
    ki = mct_aVect_indexRA(fractions,"ifrac",perrWith='quiet')
    ko = mct_aVect_indexRA(fractions,"ofrac",perrWith='quiet')
    kg = mct_aVect_indexRA(fractions,"gfrac",perrWith='quiet')
    kr = mct_aVect_indexRA(fractions,"rfrac",perrWith='quiet')
    kw = mct_aVect_indexRA(fractions,"wfrac",perrWith='quiet')
    kz = mct_aVect_indexRA(fractions,"zfrac",perrWith='quiet')
    kk = mct_aVect_indexRA(fractions,"lfrin",perrWith='quiet')

    if (ka > 0) then
       aminval = minval(fractions%rAttr(ka,:))
       amaxval = maxval(fractions%rAttr(ka,:))
    endif
    if (kl > 0) then
       lminval = minval(fractions%rAttr(kl,:))
       lmaxval = maxval(fractions%rAttr(kl,:))
    endif
    if (ko > 0) then
       ominval = minval(fractions%rAttr(ko,:))
       omaxval = maxval(fractions%rAttr(ko,:))
    endif
    if (ki > 0) then
       iminval = minval(fractions%rAttr(ki,:))
       imaxval = maxval(fractions%rAttr(ki,:))
    endif
    if (kg > 0) then
       gminval = minval(fractions%rAttr(kg,:))
       gmaxval = maxval(fractions%rAttr(kg,:))
    endif
    if (kr > 0) then
       rminval = minval(fractions%rAttr(kr,:))
       rmaxval = maxval(fractions%rAttr(kr,:))
    endif
    if (kw > 0) then
       wminval = minval(fractions%rAttr(kw,:))
       wmaxval = maxval(fractions%rAttr(kw,:))
    endif
    if (kz > 0) then
       zminval = minval(fractions%rAttr(kz,:))
       zmaxval = maxval(fractions%rAttr(kz,:))
    endif
    if (kk > 0) then
       kminval = minval(fractions%rAttr(kk,:))
       kmaxval = maxval(fractions%rAttr(kk,:))
    endif

    ncnt = 0
    maxerr = 0.0_r8
    if (kl > 0 .and. ko > 0 .and. ki > 0) then
       do n = 1,lsize
          sum = fractions%rAttr(ko,n) + fractions%rAttr(kl,n) + fractions%rAttr(ki,n)
          sminval = min(sum,sminval)
          smaxval = max(sum,smaxval)
          diff = abs(1.0_r8 - sum)
          if (diff > eps_fracsum) then
             ncnt = ncnt + 1
             maxerr = max(maxerr, diff)
             !tcx debug  write(logunit,*) trim(lstring),' err# ',ncnt, n, lsize, &
             !fractions%rAttr(ko,n),fractions%rAttr(kl,n),fractions%rAttr(ki,n),sum
          endif
       enddo
    endif

    error = .false.
    if (ncnt > 0) error = .true.
    if (aminval < 0.0_r8-eps_fracval .or. amaxval > 1.0_r8+eps_fracval) error = .true.
    if (lminval < 0.0_r8-eps_fracval .or. lmaxval > 1.0_r8+eps_fracval) error = .true.
    if (ominval < 0.0_r8-eps_fracval .or. omaxval > 1.0_r8+eps_fracval) error = .true.
    if (iminval < 0.0_r8-eps_fracval .or. imaxval > 1.0_r8+eps_fracval) error = .true.
    if (gminval < 0.0_r8-eps_fracval .or. gmaxval > 1.0_r8+eps_fracval) error = .true.
    if (rminval < 0.0_r8-eps_fracval .or. rmaxval > 1.0_r8+eps_fracval) error = .true.
    if (wminval < 0.0_r8-eps_fracval .or. wmaxval > 1.0_r8+eps_fracval) error = .true.
    if (zminval < 0.0_r8-eps_fracval .or. zmaxval > 1.0_r8+eps_fracval) error = .true.
    if (kminval < 0.0_r8-eps_fracval .or. kmaxval > 1.0_r8+eps_fracval) error = .true.

    if (error .or. seq_frac_debug > 1) then
       if (ka > 0) then
          call shr_mpi_min(aminval,tmpmin,mpicom,subname//':afrac',all=.false.)
          call shr_mpi_max(amaxval,tmpmax,mpicom,subname//':afrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' afrac min/max   = ',tmpmin,tmpmax
       endif
       if (kl > 0) then
          call shr_mpi_min(lminval,tmpmin,mpicom,subname//':lfrac',all=.false.)
          call shr_mpi_max(lmaxval,tmpmax,mpicom,subname//':lfrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' lfrac min/max   = ',tmpmin,tmpmax
       endif
       if (kg > 0) then
          call shr_mpi_min(gminval,tmpmin,mpicom,subname//':gfrac',all=.false.)
          call shr_mpi_max(gmaxval,tmpmax,mpicom,subname//':gfrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' gfrac min/max   = ',tmpmin,tmpmax
       endif
       if (ko > 0) then
          call shr_mpi_min(ominval,tmpmin,mpicom,subname//':ofrac',all=.false.)
          call shr_mpi_max(omaxval,tmpmax,mpicom,subname//':ofrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' ofrac min/max   = ',tmpmin,tmpmax
       endif
       if (ki > 0) then
          call shr_mpi_min(iminval,tmpmin,mpicom,subname//':ifrac',all=.false.)
          call shr_mpi_max(imaxval,tmpmax,mpicom,subname//':ifrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' ifrac min/max   = ',tmpmin,tmpmax
       endif
       if (kr > 0) then
          call shr_mpi_min(rminval,tmpmin,mpicom,subname//':rfrac',all=.false.)
          call shr_mpi_max(rmaxval,tmpmax,mpicom,subname//':rfrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' rfrac min/max   = ',tmpmin,tmpmax
       endif
       if (kw > 0) then
          call shr_mpi_min(wminval,tmpmin,mpicom,subname//':wfrac',all=.false.)
          call shr_mpi_max(wmaxval,tmpmax,mpicom,subname//':wfrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' wfrac min/max   = ',tmpmin,tmpmax
       endif
       if (kz > 0) then
          call shr_mpi_min(kminval,tmpmin,mpicom,subname//':zfrac',all=.false.)
          call shr_mpi_max(kmaxval,tmpmax,mpicom,subname//':zfrac',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' zfrac min/max   = ',tmpmin,tmpmax
       endif
       if (kk > 0) then
          call shr_mpi_min(kminval,tmpmin,mpicom,subname//':lfrin',all=.false.)
          call shr_mpi_max(kmaxval,tmpmax,mpicom,subname//':lfrin',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' lfrin min/max   = ',tmpmin,tmpmax
       endif
       if (kl > 0 .and. ko > 0 .and. ki > 0) then
          call shr_mpi_min(sminval,tmpmin,mpicom,subname//':sum',all=.false.)
          call shr_mpi_max(smaxval,tmpmax,mpicom,subname//':sum',all=.false.)
          if (iamroot) write(logunit,F02) trim(lstring),' sum min/max     = ',tmpmin,tmpmax
          call shr_mpi_sum(ncnt  ,tmpsum,mpicom,subname//':sum',all=.false.)
          call shr_mpi_max(maxerr,tmpmax,mpicom,subname//':sum',all=.false.)
          if (iamroot) write(logunit,F01) trim(lstring),' sum ncnt/maxerr = ',tmpsum,tmpmax
       endif
       if (error .and. .not. seq_frac_dead .and. seq_frac_abort) then
          write(logunit,F02) trim(lstring),' ERROR aborting '
          call shr_sys_abort()
       elseif (error) then
          if (iamroot) write(logunit,F02) trim(lstring),' ERROR but NOT aborting '
       endif
    endif

  end subroutine seq_frac_check

end module seq_frac_mct
