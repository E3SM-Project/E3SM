module seq_map_mod

!---------------------------------------------------------------------
!
! Purpose:
!
! General mapping routines
! including self normalizing mapping routine with optional fraction
!       
! Author: T. Craig, Jan-28-2011
!
!---------------------------------------------------------------------

  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod      ,only: CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod

  use seq_comm_mct
  use seq_mctext_mod, only : seq_mctext_gsmapextend

#ifdef USE_ESMF_LIB
  use esmf
  use esmfshr_mod
  use seq_map_esmf
#endif

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_map_init_rcfile
  public :: seq_map_init_rearrolap
  public :: seq_map_init_rearrsplit
  public :: seq_map_initvect
  public :: seq_map_map
  public :: seq_map_mapvect
  public :: seq_map_readdata
#ifdef USE_ESMF_LIB
  public :: seq_map_register
#endif

  interface seq_map_avNorm ; module procedure &
    seq_map_avNormArr, &
    seq_map_avNormAvF
  end interface

  type seq_map
    logical :: copy_only
    logical :: rearrange_only
    logical :: esmf_map
    type(mct_rearr) :: rearr
    type(mct_sMatp) :: sMatp
    !---- for cart3d
    character(CL) :: cart3d_init
    real(R8), pointer :: slon_s(:),clon_s(:),slat_s(:),clat_s(:)
    real(R8), pointer :: slon_d(:),clon_d(:),slat_d(:),clat_d(:)
    !---- for npfix
    character(CL) :: npfix_init
    integer(IN)   :: cntfound      ! number found
    integer(IN)   :: cntf_tot      ! cntfound total
    integer(IN)   :: ni,nj         ! grid size
    integer(IN)   :: mpicom        ! mpicom
    real(R8)   ,pointer :: ilon1(:)      ! lon of input grid at highest latitude
    real(R8)   ,pointer :: ilat1(:)      ! lat of input grid at highest latitude
    real(R8)   ,pointer :: ilon2(:)      ! lon of input grid at highest latitude
    real(R8)   ,pointer :: ilat2(:)      ! lat of input grid at highest latitude
    real(R8)   ,pointer :: alphafound(:) ! list of found alphas
    real(R8)   ,pointer :: betafound(:)  ! list of found betas
    integer(IN),pointer :: mfound(:)     ! list of found ms
    integer(IN),pointer :: nfound(:)     ! list of found ns
    integer(IN),pointer :: gindex(:)     ! global index 
#ifdef USE_ESMF_LIB
    !---- import and export States for this mapper object, 
    !---- routehandle is stored in the exp_state for repeated remapping use
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
#endif
  end type seq_map
  public seq_map

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  character(*),parameter :: seq_map_stroff = 'variable_unset'
  character(*),parameter :: seq_map_stron  = 'StrinG_is_ON'
  real(R8),parameter,private :: deg2rad = shr_const_pi/180.0_R8  ! deg to rads

!=======================================================================
contains
!=======================================================================

!=======================================================================
#ifdef USE_ESMF_LIB
  subroutine seq_map_register(petlist, ccsm_comp, comp, import_state, export_state)
  
    implicit none
  
    integer, pointer                  :: petlist(:)
    type(ESMF_CplComp)                :: ccsm_comp
    type(ESMF_GridComp), intent(out)  :: comp
    type(ESMF_State),    intent(out)  :: import_state, export_state
  
    integer            :: rc
  
    comp = ESMF_GridCompCreate(name="seq map comp", petList=petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create seq map comp')
    call ESMF_GridCompSetServices(comp, seq_map_esmf_register, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register atm comp')
    import_state = ESMF_StateCreate(name="seq map import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import seq map state')
    export_state = ESMF_StateCreate(name="seq map export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export seq map state')
  
  end subroutine
#endif

!=======================================================================

  subroutine seq_map_init_rearrsplit( mapper, gsmap_s, ID_s, gsmap_d, ID_d, ID_join, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout) :: mapper
    type(mct_gsmap) ,intent(in)    :: gsmap_s
    integer(IN)     ,intent(in)    :: ID_s
    type(mct_gsmap) ,intent(in)    :: gsmap_d
    integer(IN)     ,intent(in)    :: ID_d
    integer(IN)     ,intent(in)    :: ID_join
    character(len=*),intent(in),optional :: string
    !
    ! Local Variables
    !
    integer(IN) :: mpicom_s, mpicom_d, mpicom_join
    type(mct_gsmap) :: gsmap_s_join
    type(mct_gsmap) :: gsmap_d_join
    character(len=*),parameter :: subname = "(seq_map_init_rearrsplit) "

    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    mapper%copy_only = .false.
    mapper%rearrange_only = .false.
    mapper%esmf_map = .false.

    call seq_comm_setptrs(ID_join,mpicom=mpicom_join)
    mapper%mpicom = mpicom_join

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       mapper%copy_only = .true.
       if (seq_comm_iamroot(ID_join)) then
          write(logunit,'(2A,L2)') subname,' gsmaps ARE IDENTICAL, copyoption = ',mapper%copy_only
       endif
    else
       if (seq_comm_iamroot(ID_join)) write(logunit,'(2A)') subname,' gsmaps are not identical'
       mapper%rearrange_only = .true.

       call seq_comm_setptrs(ID_s ,mpicom=mpicom_s)
       call seq_comm_setptrs(ID_d ,mpicom=mpicom_d)
       call seq_comm_setptrs(ID_join,mpicom=mpicom_join)

       ! --- Extend gsmaps to join group of pes

       call seq_mctext_gsmapExtend(gsmap_s, mpicom_s, gsmap_s_join, mpicom_join, ID_join)
       call seq_mctext_gsmapExtend(gsmap_d, mpicom_d, gsmap_d_join, mpicom_join, ID_join)

       ! --- Initialize rearranger based on join gsmaps

       call seq_map_gsmapcheck(gsmap_s_join, gsmap_d_join)
       call mct_rearr_init(gsmap_s_join, gsmap_d_join, mpicom_join, mapper%rearr)

       ! --- Clean up temporary gsmaps

       call mct_gsMap_clean(gsmap_s_join)
       call mct_gsMap_clean(gsmap_d_join)

    endif

  end subroutine seq_map_init_rearrsplit

!=======================================================================

  subroutine seq_map_init_rcfile( mapper, gsmap_s, gsmap_d, mpicom, &
       maprcfile, maprcname, maprctype, samegrid, string, esmf_map)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout) :: mapper
    type(mct_gsmap) ,intent(in)    :: gsmap_s
    type(mct_gsmap) ,intent(in)    :: gsmap_d
    integer(IN)     ,intent(in)    :: mpicom
    character(len=*),intent(in)    :: maprcfile
    character(len=*),intent(in)    :: maprcname
    character(len=*),intent(in)    :: maprctype
    logical         ,intent(in)    :: samegrid
    character(len=*),intent(in),optional :: string
    logical         ,intent(in),optional :: esmf_map
    !
    ! Local Variables
    !
    integer(IN) :: ssize,dsize
    character(len=*),parameter :: subname = "(seq_map_init_rcfile) "

#ifdef USE_ESMF_LIB
    type(ESMF_GridComp)         :: comp
    type(ESMF_State)            :: imp_state, exp_state
    type(ESMF_Array)            :: gindex_s, gindex_d
    type(ESMF_Array)            :: factorArray, factorIndexArray
    integer, pointer            :: factorIndexList(:,:)
    real(ESMF_KIND_R8), pointer :: factorList(:)
    integer                     :: rc, urc, row_idx, col_idx, wgt_idx, nwgt
    logical                     :: has_weight=.true.
    integer, pointer            :: petmap(:)
#endif

    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    mapper%copy_only = .false.
    mapper%rearrange_only = .false.
    mapper%esmf_map = .false.
    mapper%mpicom = mpicom

    if (present(esmf_map)) then
       mapper%esmf_map = esmf_map
    endif

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       mapper%copy_only = .true.
    elseif (samegrid) then
       mapper%rearrange_only = .true.

       ! --- Initialize rearranger
       call seq_map_gsmapcheck(gsmap_s, gsmap_d)
       call mct_rearr_init(gsmap_s, gsmap_d, mpicom, mapper%rearr)

    else

       ! --- Initialize Smatp
       call shr_mct_sMatPInitnc(mapper%sMatp,gsMap_s,gsMap_d, &
            trim(maprcfile),trim(maprcname),trim(maprctype),mpicom)

       if (mapper%esmf_map) then
#ifdef USE_ESMF_LIB
          !------------------------------------------------------
          ! ESMF:
          ! Set up routehandle
          !   This section of the code handles the creation of routehandle
          ! for ESMF SMM communication calls. First, the existing MCT
          ! sparsematmul descriptor is reused (although this part can be
          ! done completely without MCT, doing a distributed IO from maprcfile).
          ! The weight matrix and indicies vectors are wrapped into ESMF Arrays
          ! based on simple block distributed DistGrid. 
          !   These Arrays are passed into a coupler component set up during
          ! CESM driver initialization and stored in global seq_comm_type array
          ! at the driver level. The coupler component is indexed by CPLID in
          ! the seq_comm_type array. In the coupler component initialization phase
          ! the SMM routehandle is computed and attached to export state of
          ! a specific mapper object. Routehandle in a mapper object can be reused
          ! for the same regridding model component pair.
          !
          ! Fei Liu
          !------------------------------------------------------
          ! has_weight can be controlled through namelist to determine
          ! if weights are computed offline or to be computed online
          ! during initialization.
          if(has_weight) then
            ! Retrieve weights and indicies from mapper object
            row_idx = mct_sMat_indexIA(mapper%sMatp%matrix, 'grow')
            col_idx = mct_sMat_indexIA(mapper%sMatp%matrix, 'gcol')
            wgt_idx = mct_sMat_indexRA(mapper%sMatp%matrix, 'weight')
            nwgt = size(mapper%sMatp%matrix%data%rattr(wgt_idx, :))
            !write(logunit,*) trim(string), row_idx, col_idx, wgt_idx, nwgt
            allocate(factorIndexList(2, nwgt), factorList(nwgt), stat=rc)
            if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
            factorIndexList(1,:) = mapper%sMatp%matrix%data%iattr(col_idx, :)
            factorIndexList(2,:) = mapper%sMatp%matrix%data%iattr(row_idx, :)
            factorList = mapper%sMatp%matrix%data%rattr(wgt_idx,:)

            ! Get coupler pet map in global setting
            call seq_comm_petlist(CPLID, petmap)

            ! Set up temporary arrays to compute ESMF SMM routes.
            gindex_s = mct2esmf_create(gsMap_s, mpicom=mpicom, petmap=petmap, &
                name="gindex_s", rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            gindex_d = mct2esmf_create(gsMap_d, mpicom=mpicom, petmap=petmap, &
                name="gindex_d", rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            factorArray = mct2esmf_create(factorList, mpicom=mpicom, petmap=petmap, &
                name="factorArray", rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            factorIndexArray = mct2esmf_create(factorIndexList, mpicom=mpicom, &
                petmap=petmap, name="factorIndexArray", rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            ! Get the coupler component
            ! Create mapper specific import and export States
            call seq_comm_getcompstates(CPLID, comp=comp)

            mapper%imp_state = ESMF_StateCreate(name="import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
            if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import atm state')
            mapper%exp_state = ESMF_StateCreate(name="export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
            if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export atm state')

            ! Attach Arrays to the States
            call ESMF_StateAdd(mapper%imp_state, (/gindex_s/), rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            call ESMF_StateAdd(mapper%exp_state, (/gindex_d/), rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            call ESMF_StateAdd(mapper%exp_state, (/factorArray/), rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            call ESMF_StateAdd(mapper%exp_state, (/factorIndexArray/), rc=rc)
            if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            ! Call into ESMF init method
            call ESMF_GridCompInitialize(comp, importState=mapper%imp_state, exportState=mapper%exp_state, &
                 userRc=urc, rc=rc)
            if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
            if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

            deallocate(factorIndexList, factorList)

          else
            ! Compute the routehandle online here:
            ! 1. Read in source and destination Grid description files
            ! 2. Construct Grids or Meshes and dummy Fields on top
            ! 3. Call FieldRegridStore to compute routehandle
            ! 4. Save routehandle in mapper
          endif
#else
          call shr_sys_abort(subname//' ERROR: esmf SMM not allowed without USE_ESMF_LIB')
#endif
       endif  ! esmf_map

    endif

  end subroutine seq_map_init_rcfile

!=======================================================================

  subroutine seq_map_init_rearrolap( mapper, gsmap_s, gsmap_d, mpicom, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout) :: mapper
    type(mct_gsmap) ,intent(in)    :: gsmap_s
    type(mct_gsmap) ,intent(in)    :: gsmap_d
    integer(IN)     ,intent(in)    :: mpicom
    character(len=*),intent(in),optional :: string
    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(seq_map_init_rearrolap) "

    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    mapper%copy_only = .false.
    mapper%rearrange_only = .false.
    mapper%mpicom = mpicom

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       mapper%copy_only = .true.
    else
       mapper%rearrange_only = .true.

       ! --- Initialize rearranger
       call seq_map_gsmapcheck(gsmap_s, gsmap_d)
       call mct_rearr_init(gsmap_s, gsmap_d, mpicom, mapper%rearr)

    endif

  end subroutine seq_map_init_rearrolap

!=======================================================================

  subroutine seq_map_map( mapper, av_s, av_d, fldlist, norm, avwts_s, avwtsfld_s, &
                          string, msgtag )

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout):: mapper
    type(mct_aVect) ,intent(in)   :: av_s
    type(mct_aVect) ,intent(inout):: av_d
    character(len=*),intent(in),optional :: fldlist
    logical         ,intent(in),optional :: norm
    type(mct_aVect) ,intent(in),optional :: avwts_s
    character(len=*),intent(in),optional :: avwtsfld_s
    character(len=*),intent(in),optional :: string
    integer(IN)     ,intent(in),optional :: msgtag
    !
    ! Local Variables
    !
    logical :: lnorm
    integer(IN),save :: ltag    ! message tag for rearrange
    character(len=*),parameter :: subname = "(seq_map_map) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (present(msgtag)) then
       ltag = msgtag
    else
       ltag = 2000
    endif

    if (present(avwts_s) .and. .not. present(avwtsfld_s)) then
       write(logunit,*) subname,' ERROR: avwts_s present but avwtsfld_s not'
       call shr_sys_abort(subname//' ERROR: avwts present')
    endif
    if (.not. present(avwts_s) .and. present(avwtsfld_s)) then
       write(logunit,*) subname,' ERROR: avwtsfld_s present but avwts_s not'
       call shr_sys_abort(subname//' ERROR: avwtsfld present')
    endif

    if (mapper%copy_only) then
       !-------------------------------------------
       ! COPY data
       !-------------------------------------------
       if (present(fldlist)) then
          call mct_aVect_copy(aVin=av_s,aVout=av_d,rList=fldlist,vector=mct_usevector)
       else
          call mct_aVect_copy(aVin=av_s,aVout=av_d,vector=mct_usevector)
       endif

    else if (mapper%rearrange_only) then
       !-------------------------------------------
       ! REARRANGE data
       !-------------------------------------------
       if (present(fldlist)) then
          call mct_rearr_rearrange_fldlist(av_s, av_d, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall, fldlist=fldlist)
       else
          call mct_rearr_rearrange(av_s, av_d, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall)
       endif

    else
       !-------------------------------------------
       ! MAP data
       !-------------------------------------------
       if (present(avwts_s)) then
          if (present(fldlist)) then
             call seq_map_avNorm(mapper, av_s, av_d, avwts_s, trim(avwtsfld_s), &
                  rList=fldlist, norm=lnorm)
          else
             call seq_map_avNorm(mapper, av_s, av_d, avwts_s, trim(avwtsfld_s), &
                  norm=lnorm)
          endif
       else
          if (present(fldlist)) then
             call seq_map_avNorm(mapper, av_s, av_d, rList=fldlist, norm=lnorm)
          else
             call seq_map_avNorm(mapper, av_s, av_d, norm=lnorm)
          endif
       endif
    end if

  end subroutine seq_map_map

!=======================================================================
  subroutine seq_map_initvect( mapper, type, dom_s, dom_d, gsmap_s, ni, nj, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout):: mapper
    character(len=*),intent(in)   :: type
    type(mct_gGrid) ,intent(in)   :: dom_s
    type(mct_gGrid) ,intent(inout):: dom_d
    type(mct_gsMap) ,intent(in),optional :: gsmap_s
    integer(IN)     ,intent(in),optional :: ni
    integer(IN)     ,intent(in),optional :: nj
    character(len=*),intent(in),optional :: string
    !
    ! Local Variables
    !
    integer(IN) :: klon, klat, lsize, n
    logical :: lnorm
    character(len=CL) :: lstring
    character(len=*),parameter :: subname = "(seq_map_initvect) "

    lstring = ' '
    if (present(string)) then
!       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
       lstring = trim(string)
    endif

    if (trim(type) == 'npfix') then

       if (mapper%npfix_init == trim(seq_map_stron)) return

       if (.not.present(gsmap_s)) then
          call shr_sys_abort(trim(subname)//' ERROR gsmap_s required for npfix')
       endif
       if (.not.present(ni)) then
          call shr_sys_abort(trim(subname)//' ERROR ni required for npfix')
       endif
       if (.not.present(nj)) then
          call shr_sys_abort(trim(subname)//' ERROR nj required for npfix')
       endif
       mapper%ni = ni
       mapper%nj = nj
       call map_npFixNew4R(mapper,'init',gsmapi=gsmap_s,domi=dom_s,domo=dom_d)
       mapper%npfix_init = trim(seq_map_stron)

    elseif (trim(type(1:6)) == 'cart3d') then
       if (mapper%cart3d_init == trim(seq_map_stron)) return

       !--- compute these up front for vector mapping ---
       lsize = mct_aVect_lsize(dom_s%data)
       allocate(mapper%slon_s(lsize),mapper%clon_s(lsize), &
                mapper%slat_s(lsize),mapper%clat_s(lsize))
       klon = mct_aVect_indexRa(dom_s%data, "lon" )
       klat = mct_aVect_indexRa(dom_s%data, "lat" )
       do n = 1,lsize
          mapper%slon_s(n) = sin(dom_s%data%rAttr(klon,n)*deg2rad)
          mapper%clon_s(n) = cos(dom_s%data%rAttr(klon,n)*deg2rad)
          mapper%slat_s(n) = sin(dom_s%data%rAttr(klat,n)*deg2rad)
          mapper%clat_s(n) = cos(dom_s%data%rAttr(klat,n)*deg2rad)
       enddo
       
       lsize = mct_aVect_lsize(dom_d%data)
       allocate(mapper%slon_d(lsize),mapper%clon_d(lsize), &
                mapper%slat_d(lsize),mapper%clat_d(lsize))
       klon = mct_aVect_indexRa(dom_d%data, "lon" )
       klat = mct_aVect_indexRa(dom_d%data, "lat" )
       do n = 1,lsize
          mapper%slon_d(n) = sin(dom_d%data%rAttr(klon,n)*deg2rad)
          mapper%clon_d(n) = cos(dom_d%data%rAttr(klon,n)*deg2rad)
          mapper%slat_d(n) = sin(dom_d%data%rAttr(klat,n)*deg2rad)
          mapper%clat_d(n) = cos(dom_d%data%rAttr(klat,n)*deg2rad)
       enddo
       mapper%cart3d_init = trim(seq_map_stron)
    endif

  end subroutine seq_map_initvect

!=======================================================================
  subroutine seq_map_mapvect( mapper, type, av_s, av_d, fldu, fldv, norm, avwts_s, avwtsfld_s, string )

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout):: mapper
    character(len=*),intent(in)   :: type
    type(mct_aVect) ,intent(in)   :: av_s
    type(mct_aVect) ,intent(inout):: av_d
    character(len=*),intent(in)   :: fldu
    character(len=*),intent(in)   :: fldv
    logical         ,intent(in),optional :: norm
    type(mct_aVect) ,intent(in),optional :: avwts_s
    character(len=*),intent(in),optional :: avwtsfld_s
    character(len=*),intent(in),optional :: string
    !
    ! Local Variables
    !
    logical :: lnorm
    character(len=CL) :: lstring
    character(len=*),parameter :: subname = "(seq_map_mapvect) "
    !-----------------------------------------------------

    lstring = ' '
    if (present(string)) then
!       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
       lstring = trim(string)
    endif

    if (mapper%copy_only .or. mapper%rearrange_only) then
       return
    endif

    if ((present(avwts_s) .and. .not. present(avwtsfld_s)) .or. &
        (present(avwtsfld_s) .and. .not. present(avwts_s))) then
       call shr_sys_abort(trim(subname)//' ERROR avwts_s avwtsfld_s consistency')
    endif

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (trim(type) == 'npfix') then
       if (mapper%npfix_init /= trim(seq_map_stron)) then
          call shr_sys_abort(trim(subname)//' ERROR: npfix not initialized '//trim(lstring))
       endif
       call map_npFixNew4R(mapper,'run',buni=av_s,buno=av_d,fld1=fldu,fld2=fldv)
    elseif (trim(type(1:6)) == 'cart3d') then
       if (mapper%cart3d_init /= trim(seq_map_stron)) then
          call shr_sys_abort(trim(subname)//' ERROR: cart3d not initialized '//trim(lstring))
       endif
       if (present(avwts_s) .and. present(avwtsfld_s)) then
          call seq_map_cart3d(mapper, type, av_s, av_d, fldu, fldv, norm=lnorm, &
               avwts_s=avwts_s, avwtsfld_s=avwtsfld_s, string=string)
       else
          call seq_map_cart3d(mapper, type, av_s, av_d, fldu, fldv, norm=lnorm, string=string)
       endif
    elseif (trim(type) == 'none') then
       call seq_map_map(mapper, av_s, av_d, fldlist=trim(fldu)//':'//trim(fldv), norm=lnorm)
    else
       write(logunit,*) subname,' ERROR: type unsupported ',trim(type)
       call shr_sys_abort(trim(subname)//' ERROR in type='//trim(type))
    end if

  end subroutine seq_map_mapvect

!=======================================================================
  subroutine seq_map_cart3d( mapper, type, av_s, av_d, fldu, fldv, norm, avwts_s, avwtsfld_s, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout):: mapper
    character(len=*),intent(in)   :: type
    type(mct_aVect) ,intent(in)   :: av_s
    type(mct_aVect) ,intent(inout):: av_d
    character(len=*),intent(in)   :: fldu
    character(len=*),intent(in)   :: fldv
    logical         ,intent(in),optional :: norm
    type(mct_aVect) ,intent(in),optional :: avwts_s
    character(len=*),intent(in),optional :: avwtsfld_s
    character(len=*),intent(in),optional :: string
    !
    ! Local Variables
    !
    integer           :: lsize
    logical           :: lnorm
    integer           :: ku,kv,kux,kuy,kuz,n
    real(r8)          :: ue,un,ur,ux,uy,uz,speed
    real(r8)          :: urmaxl,urmax,uravgl,uravg,spavgl,spavg
    type(mct_aVect)   :: av3_s, av3_d
    integer(in)       :: mpicom,my_task,ierr,urcnt,urcntl
    character(len=*),parameter :: subname = "(seq_map_cart3d) "

    lnorm = .true.
    if (present(norm)) then
       lnorm=norm
    endif

    mpicom = mapper%mpicom

    if ((present(avwts_s) .and. .not. present(avwtsfld_s)) .or. &
        (present(avwtsfld_s) .and. .not. present(avwts_s))) then
       call shr_sys_abort(trim(subname)//' ERROR avwts_s avwtsfld_s consistency')
    endif

    ku = mct_aVect_indexRA(av_s, trim(fldu), perrwith='quiet')
    kv = mct_aVect_indexRA(av_s, trim(fldv), perrwith='quiet')

    if (ku /= 0 .and. kv /= 0) then
       lsize = mct_aVect_lsize(av_s)
       call mct_avect_init(av3_s,rList='ux:uy:uz',lsize=lsize)

       lsize = mct_aVect_lsize(av_d)
       call mct_avect_init(av3_d,rList='ux:uy:uz',lsize=lsize)

       kux = mct_aVect_indexRA(av3_s,'ux')
       kuy = mct_aVect_indexRA(av3_s,'uy')
       kuz = mct_aVect_indexRA(av3_s,'uz')
       lsize = mct_aVect_lsize(av_s)
       do n = 1,lsize
          ur = 0.0_r8
          ue = av_s%rAttr(ku,n)
          un = av_s%rAttr(kv,n)
          ux = mapper%clon_s(n)*mapper%clat_s(n)*ur - &
               mapper%clon_s(n)*mapper%slat_s(n)*un - &
               mapper%slon_s(n)*ue
          uy = mapper%slon_s(n)*mapper%clon_s(n)*ur - &
               mapper%slon_s(n)*mapper%slat_s(n)*un + &
               mapper%clon_s(n)*ue
          uz = mapper%slat_s(n)*ur + &
               mapper%clat_s(n)*un
          av3_s%rAttr(kux,n) = ux
          av3_s%rAttr(kuy,n) = uy
          av3_s%rAttr(kuz,n) = uz
       enddo

       if (present(avwts_s) .and. present(avwtsfld_s)) then
          call seq_map_map(mapper, av3_s, av3_d, norm=lnorm, avwts_s=avwts_s, avwtsfld_s=avwtsfld_s)
       else
          call seq_map_map(mapper, av3_s, av3_d, norm=lnorm)
       endif

       kux = mct_aVect_indexRA(av3_d,'ux')
       kuy = mct_aVect_indexRA(av3_d,'uy')
       kuz = mct_aVect_indexRA(av3_d,'uz')
       lsize = mct_aVect_lsize(av_d)
       urmaxl = -1.0_r8
       uravgl = 0.0_r8
       urcntl = 0
       spavgl = 0.0_r8
       do n = 1,lsize
          ux = av3_d%rAttr(kux,n)
          uy = av3_d%rAttr(kuy,n)
          uz = av3_d%rAttr(kuz,n)
          ue = -mapper%slon_d(n)          *ux + &
                mapper%clon_d(n)          *uy
          un = -mapper%clon_d(n)*mapper%slat_d(n)*ux - &
                mapper%slon_d(n)*mapper%slat_d(n)*uy + &
                mapper%clat_d(n)*uz
          ur =  mapper%clon_d(n)*mapper%clat_d(n)*ux + &
                mapper%slon_d(n)*mapper%clat_d(n)*uy - &
                mapper%slat_d(n)*uz
          speed = sqrt(ur*ur + ue*ue + un*un)
          if (trim(type) == 'cart3d_diag' .or. trim(type) == 'cart3d_uvw_diag') then
             if (speed /= 0.0_r8) then
                urmaxl = max(urmaxl,abs(ur))
                uravgl = uravgl + abs(ur)
                spavgl = spavgl + speed
                urcntl = urcntl + 1
             endif
          endif
          if (type(1:10) == 'cart3d_uvw') then
             !--- this adds ur to ue and un, while preserving u/v angle and total speed ---
             if (un == 0.0_R8) then
                !--- if ue is also 0.0 then just give speed to ue, this is arbitrary ---
                av_d%rAttr(ku,n) = sign(speed,ue)
                av_d%rAttr(kv,n) = 0.0_r8
             else if (ue == 0.0_R8) then
                av_d%rAttr(ku,n) = 0.0_r8
                av_d%rAttr(kv,n) = sign(speed,un)
             else
                av_d%rAttr(ku,n) = sign(speed/sqrt(1.0_r8 + ((un*un)/(ue*ue))),ue)
                av_d%rAttr(kv,n) = sign(speed/sqrt(1.0_r8 + ((ue*ue)/(un*un))),un)
             endif
          else
             !--- this ignores ur ---
             av_d%rAttr(ku,n) = ue
             av_d%rAttr(kv,n) = un
          endif
       enddo
       if (trim(type) == 'cart3d_diag' .or. trim(type) == 'cart3d_uvw_diag') then
          call mpi_comm_rank(mpicom,my_task,ierr)
          call shr_mpi_max(urmaxl,urmax,mpicom,'urmax')
          call shr_mpi_sum(uravgl,uravg,mpicom,'uravg')
          call shr_mpi_sum(spavgl,spavg,mpicom,'spavg')
          call shr_mpi_sum(urcntl,urcnt,mpicom,'urcnt')
          if (my_task == 0 .and. urcnt > 0) then
             uravg = uravg / urcnt
             spavg = spavg / urcnt
             write(logunit,*) trim(subname),' cart3d uravg,urmax,spavg = ',uravg,urmax,spavg
          endif
       endif

      call mct_avect_clean(av3_s)
      call mct_avect_clean(av3_d)

   endif  ! ku,kv

  end subroutine seq_map_cart3d
!=======================================================================

  subroutine seq_map_readdata(maprcfile, maprcname, mpicom, ID, &
         ni_s, nj_s, av_s, gsmap_s, avfld_s, filefld_s, &
         ni_d, nj_d, av_d, gsmap_d, avfld_d, filefld_d, string)

    !--- lifted from work by J Edwards, April 2011

    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    use pio, only : pio_openfile, pio_closefile, pio_read_darray, pio_inq_dimid, &
       pio_inq_dimlen, pio_inq_varid, file_desc_t, io_desc_t, iosystem_desc_t, &
       var_desc_t, pio_int, pio_get_var, pio_double, pio_initdecomp, pio_freedecomp
    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    character(len=*),intent(in)    :: maprcfile
    character(len=*),intent(in)    :: maprcname
    integer(IN)     ,intent(in)    :: mpicom
    integer(IN)     ,intent(in)    :: ID
    integer(IN)     ,intent(out)  ,optional :: ni_s
    integer(IN)     ,intent(out)  ,optional :: nj_s
    type(mct_avect) ,intent(inout),optional :: av_s
    type(mct_gsmap) ,intent(in)   ,optional :: gsmap_s
    character(len=*),intent(in)   ,optional :: avfld_s
    character(len=*),intent(in)   ,optional :: filefld_s
    integer(IN)     ,intent(out)  ,optional :: ni_d
    integer(IN)     ,intent(out)  ,optional :: nj_d
    type(mct_avect) ,intent(inout),optional :: av_d
    type(mct_gsmap) ,intent(in)   ,optional :: gsmap_d
    character(len=*),intent(in)   ,optional :: avfld_d
    character(len=*),intent(in)   ,optional :: filefld_d
    character(len=*),intent(in)   ,optional :: string
    !
    ! Local Variables
    !
    type(iosystem_desc_t), pointer :: pio_subsystem
    integer(IN)       :: pio_iotype
    type(file_desc_t) :: File    ! PIO file pointer
    type(io_desc_t)   :: iodesc  ! PIO parallel io descriptor
    integer(IN)       :: rcode   ! pio routine return code
    type(var_desc_t)  :: vid     ! pio variable  ID
    integer(IN)       :: did     ! pio dimension ID
    integer(IN)       :: na      ! size of source domain
    integer(IN)       :: nb      ! size of destination domain
    integer(IN)       :: i       ! index
    integer(IN)       :: mytask  ! my task
    integer(IN), pointer :: dof(:)    ! DOF pointers for parallel read
    character(len=256):: fileName
    character(len=64) :: lfld_s, lfld_d, lfile_s, lfile_d
    character(*),parameter :: areaAV_field = 'aream'
    character(*),parameter :: areafile_s   = 'area_a'
    character(*),parameter :: areafile_d   = 'area_b'
    character(len=*),parameter :: subname  = "(seq_map_readdata) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
!       call shr_sys_flush(logunit)
    endif

    call MPI_COMM_RANK(mpicom,mytask,rcode)

    lfld_s = trim(areaAV_field)
    if (present(avfld_s)) then
       lfld_s = trim(avfld_s)
    endif

    lfld_d = trim(areaAV_field)
    if (present(avfld_d)) then
       lfld_s = trim(avfld_d)
    endif

    lfile_s = trim(areafile_s)
    if (present(filefld_s)) then
       lfile_s = trim(filefld_s)
    endif

    lfile_d = trim(areafile_d)
    if (present(filefld_d)) then
       lfile_d = trim(filefld_d)
    endif

    call I90_allLoadF(trim(maprcfile),0,mpicom,rcode)
    if(rcode /= 0) then
       write(logunit,*)"Cant find maprcfile file ",trim(maprcfile)
       call shr_sys_abort(trim(subname)//"i90_allLoadF File Not Found")
    endif

    call i90_label(trim(maprcname),rcode)
    if(rcode /= 0) then
       write(logunit,*)"Cant find label ",maprcname
       call shr_sys_abort(trim(subname)//"i90_label Not Found")
    endif

    call i90_gtoken(filename,rcode)
    if(rcode /= 0) then
       write(logunit,*)"Error reading token ",filename
       call shr_sys_abort(trim(subname)//"i90_gtoken Error on filename read")
    endif

    pio_subsystem => shr_pio_getiosys(ID)
    pio_iotype = shr_pio_getiotype(ID)

    rcode = pio_openfile(pio_subsystem, File, pio_iotype, filename)

    if (present(ni_s)) then 
       rcode = pio_inq_dimid (File, 'ni_a', did)  ! number of lons in input grid
       rcode = pio_inq_dimlen(File, did  , ni_s)
    end if
    if(present(nj_s)) then
       rcode = pio_inq_dimid (File, 'nj_a', did)  ! number of lats in input grid
       rcode = pio_inq_dimlen(File, did  , nj_s)
    end if
    if(present(ni_d)) then
       rcode = pio_inq_dimid (File, 'ni_b', did)  ! number of lons in output grid
       rcode = pio_inq_dimlen(File, did  , ni_d)
    end if
    if(present(nj_d)) then
       rcode = pio_inq_dimid (File, 'nj_b', did)  ! number of lats in output grid
       rcode = pio_inq_dimlen(File, did  , nj_d)
    endif

    !--- read and load area_a ---
    if (present(av_s)) then
       if (.not.present(gsmap_s)) then
          call shr_sys_abort(trim(subname)//' ERROR av_s must have gsmap_s')
       endif
       rcode = pio_inq_dimid (File, 'n_a', did)  ! size of  input vector
       rcode = pio_inq_dimlen(File, did  , na)
       i = mct_avect_indexra(av_s, trim(lfld_s))
       call mct_gsmap_OrderedPoints(gsMap_s, mytask, dof)
       call pio_initdecomp(pio_subsystem, pio_double, (/na/), dof, iodesc)
       deallocate(dof)
       rcode = pio_inq_varid(File,trim(lfile_s),vid)
       call pio_read_darray(File, vid, iodesc, av_s%rattr(i,:), rcode)
       call pio_freedecomp(File,iodesc)
    end if

    !--- read and load area_b ---
    if (present(av_d)) then
       if (.not.present(gsmap_d)) then
          call shr_sys_abort(trim(subname)//' ERROR av_d must have gsmap_d')
       endif
       rcode = pio_inq_dimid (File, 'n_b', did)  ! size of output vector
       rcode = pio_inq_dimlen(File, did  , nb)
       i = mct_avect_indexra(av_d, trim(lfld_d))
       call mct_gsmap_OrderedPoints(gsMap_d, mytask, dof)
       call pio_initdecomp(pio_subsystem, pio_double, (/nb/), dof, iodesc)
       deallocate(dof)
       rcode = pio_inq_varid(File,trim(lfile_d),vid)
       call pio_read_darray(File, vid, iodesc, av_d%rattr(i,:), rcode)
       call pio_freedecomp(File,iodesc)
    endif


    call pio_closefile(File)

  end subroutine seq_map_readdata

!=======================================================================

  subroutine seq_map_avNormAvF(mapper, av_i, av_o, avf_i, avfifld, rList, norm)

    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_map)   , intent(inout) :: mapper! mapper
    type(mct_aVect) , intent(in)    :: av_i  ! input 
    type(mct_aVect) , intent(inout) :: av_o  ! output
    type(mct_aVect) , intent(in)    :: avf_i  ! extra src "weight"
    character(len=*), intent(in)    :: avfifld ! field name in avf_i
    character(len=*), intent(in),optional :: rList ! fields list
    logical         , intent(in),optional :: norm  ! normalize at end
    !
    integer(IN) :: lsize_i, lsize_f, lsize_o, kf, j
    real(r8),allocatable :: frac_i(:),frac_o(:)
    logical :: lnorm
    character(*),parameter :: subName = '(seq_map_avNormAvF) '
    !-----------------------------------------------------

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    lsize_i = mct_aVect_lsize(av_i)
    lsize_f = mct_aVect_lsize(avf_i)

    if (lsize_i /= lsize_f) then
       write(logunit,*) subname,' ERROR: lsize_i ne lsize_f ',lsize_i,lsize_f
       call shr_sys_abort(subname//' ERROR size_i ne lsize_f')
    endif

    !--- extract frac_i field from avf_i to pass to seq_map_avNormArr ---
    allocate(frac_i(lsize_i))
    do j = 1,lsize_i
       kf = mct_aVect_indexRA(avf_i,trim(avfifld))
       frac_i(j) = avf_i%rAttr(kf,j)
    enddo

    if (present(rList)) then
       call seq_map_avNormArr(mapper, av_i, av_o, frac_i, rList=rList, norm=lnorm)
    else
       call seq_map_avNormArr(mapper, av_i, av_o, frac_i, norm=lnorm)
    endif

    deallocate(frac_i)

  end subroutine seq_map_avNormAvF

!=======================================================================

  subroutine seq_map_avNormArr(mapper, av_i, av_o, norm_i, rList, norm)

    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_map)   , intent(inout) :: mapper! mapper
    type(mct_aVect) , intent(in)    :: av_i  ! input 
    type(mct_aVect) , intent(inout) :: av_o  ! output
    real(r8)        , intent(in), optional :: norm_i(:)  ! source "weight"
    character(len=*), intent(in), optional :: rList ! fields list
    logical         , intent(in), optional :: norm  ! normalize at end
    !
    ! Local variables
    !
    type(mct_sMatp)        :: sMatp ! sMat
    type(mct_aVect)        :: avp_i , avp_o
    integer(IN)            :: i,j,ier,kf
    integer(IN)            :: lsize_i,lsize_o
    real(r8)               :: normval
    character(CX)          :: lrList
    logical                :: lnorm
    character(*),parameter :: subName = '(seq_map_avNormArr) '
    character(len=*),parameter :: ffld = 'norm8wt'  ! want something unique
#ifdef USE_ESMF_LIB
    type(ESMF_Array)    :: array_s, array_d, temparray_s, temparray_d
    type(ESMF_DistGrid) :: distgrid
    type(ESMF_GridComp) :: comp
    integer             :: rc, urc
#endif
    !-----------------------------------------------------

    sMatp   = mapper%sMatp
    lsize_i = mct_aVect_lsize(av_i)
    lsize_o = mct_aVect_lsize(av_o)

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (present(norm_i) .and..not.lnorm) then
       write(logunit,*) subname,' ERROR: norm_i and norm = false'
       call shr_sys_abort(subname//' ERROR norm_i and norm = false')
    endif

    if (present(norm_i)) then
       if (size(norm_i) /= lsize_i) then
          write(logunit,*) subname,' ERROR: size(norm_i) ne lsize_i ',size(norm_i),lsize_i
          call shr_sys_abort(subname//' ERROR size(norm_i) ne lsize_i')
       endif
    endif

    !--- create temporary avs for mapping ---

    if (present(rList)) then
       call mct_aVect_init(avp_i, rList=trim( rList)//':'//ffld, lsize=lsize_i)
       call mct_aVect_init(avp_o, rList=trim( rList)//':'//ffld, lsize=lsize_o)
    else
       lrList = trim(mct_aVect_exportRList2c(av_i))
       call mct_aVect_init(avp_i, rList=trim(lrList)//':'//ffld, lsize=lsize_i)
       lrList = trim(mct_aVect_exportRList2c(av_o))
       call mct_aVect_init(avp_o, rList=trim(lrList)//':'//ffld, lsize=lsize_o)
    endif

    !--- copy av_i to avp_i and set ffld value to 1.0
    !--- then multiply all fields by norm_i if norm_i exists 
    !--- this will do the right thing for the norm_i normalization 

    call mct_aVect_copy(aVin=av_i, aVout=avp_i, VECTOR=mct_usevector)
    kf = mct_aVect_indexRA(avp_i,ffld)
    do j = 1,lsize_i
       avp_i%rAttr(kf,j) = 1.0_r8
    enddo

    if (present(norm_i)) then
       do j = 1,lsize_i
          avp_i%rAttr(:,j) = avp_i%rAttr(:,j)*norm_i(j)
       enddo
    endif

    !--- map ---

    if (mapper%esmf_map) then
#ifdef USE_ESMF_LIB
       !--- ESMF SMM ---
       !   The 2nd part of the code in seq_map_mod.F90 handles the execution of ESMF SMM. Based
       ! on the number of input and output physical fields, temporary input and output ESMF Arrays
       ! are created and reference data pointer directly in input/output MCT AttributeVectors. These
       ! Arrays are attached to the import and export State retrieved from the mapper object of
       ! the regridding model component pair prior to the component run method. 
       ! Inside the component run method, a SMM execution is performed with the input/output Arrays
       ! and the routehandle. Afterwards, the temporary input/output Array are removed from 
       ! the States and destroyed.

       ! Create the temporary input array and attach to import State
       !
       ! Query DistGrid
       call ESMF_StateGet(mapper%imp_state, itemName='temparray_s', array=temparray_s, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_ArrayGet(temparray_s, distgrid=distgrid, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Create the temporary input array
       array_s = ESMF_ArrayCreate(distgrid=distgrid, farrayPtr=avp_i%rattr, &
           distgridToArrayMap=(/2/), &
           name='array_s', rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_AttributeSet(array_s, name='mct_names', &
           value=trim(mct_aVect_exportRList2c(avp_i)), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Attach to import State
       call ESMF_StateAdd(mapper%imp_state, (/array_s/), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Create the temporary output array and attach to export State
       call ESMF_StateGet(mapper%exp_state, itemName='temparray_d', array=temparray_d, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_ArrayGet(temparray_d, distgrid=distgrid, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       array_d = ESMF_ArrayCreate(distgrid=distgrid, farrayPtr=avp_o%rattr, &
           distgridToArrayMap=(/2/), &
           name='array_d', rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_AttributeSet(array_d, name='mct_names', &
           value=trim(mct_aVect_exportRList2c(avp_o)), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(mapper%exp_state, (/array_d/), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Query coupler gridcomp component and call into ESMF comp run method for SMM execution
       call seq_comm_getcompstates(CPLID, comp=comp)
       call ESMF_GridCompRun(comp, importState=mapper%imp_state, exportState=mapper%exp_state, &
           userRc=urc, rc=rc)
       if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Remove the input/output arrays from import/export States and destroy them
       call ESMF_StateRemove(mapper%imp_state, itemName='array_s', rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_ArrayDestroy(array_s, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateRemove(mapper%exp_state, itemName='array_d', rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_ArrayDestroy(array_d, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
#else
       call shr_sys_abort(subname//' ERROR: esmf SMM not allowed without USE_ESMF_LIB')
#endif

    else
       ! MCT based SMM
       call mct_sMat_avMult(avp_i, sMatp, avp_o, VECTOR=mct_usevector)
    endif


    !--- renormalize avp_o by mapped norm_i  ---

    if (lnorm) then
       do j = 1,lsize_o
          kf = mct_aVect_indexRA(avp_o,ffld)
          normval = avp_o%rAttr(kf,j)
          if (normval /= 0.0_r8) then
             normval = 1.0_r8/normval
          endif
          avp_o%rAttr(:,j) = avp_o%rAttr(:,j)*normval
       enddo
    endif

    !--- copy back into av_o and we are done ---

    call mct_aVect_copy(aVin=avp_o, aVout=av_o, VECTOR=mct_usevector)

    call mct_aVect_clean(avp_i)
    call mct_aVect_clean(avp_o)

 end subroutine seq_map_avNormArr

!=======================================================================

 subroutine map_npFixNew4R(mapper,mode,buni,buno,fld1,fld2,gsmapi,domi,domo)

   !===============================================================================
   !    Correct the north pole mapping of velocity fields from the atm to ocn
   !    grids.  This assumes the input grid is a regular lat/lon with the north
   !    pole surrounded by the last latitude line in the input array.  The
   !    longitudes in the last latitude must be ordered and equally spaced.
   !
   !    4R is a low memory version of 3R.
   !    This version (New4R) is the same as 3R except it uses a lot less memory
   !    and is a bit faster.  Like 3R, it saves data between calls and so
   !    assumes the input grid remains constant for all calls.  This is bfb
   !    with 3R on bluevista as of 2/15/07.
   !
   !    !REVISION HISTORY:
   !    2007-Feb-12 - T. Craig -- modified New3R to reduce memory
   !    2007-Apr-27 - M. Vertenstein - implemented in sequential system
   !===============================================================================
   
#include <mpif.h>  
   ! 
   ! Arguments
   !
   type(seq_map)  ,intent(inout):: mapper  ! mapper
   character(*)   ,intent(in)   :: mode    ! init or run
   type(mct_Avect),intent(in)   ,optional:: buni    ! input  attribute vec
   type(mct_Avect),intent(inout),optional:: buno    ! output attribute vec
   character(*)   ,intent(in)   ,optional:: fld1    ! name of first input field
   character(*)   ,intent(in)   ,optional:: fld2    ! name of second input field
   type(mct_gsMap),intent(in)   ,optional:: gsmapi  ! input gsmap
   type(mct_gGrid),intent(in)   ,optional:: domi    ! input domain 
   type(mct_gGrid),intent(in)   ,optional:: domo    ! output domain 
   !
   ! Local Variables
   !
   integer(IN)  :: n,m                     ! generic indices
   integer(IN)  :: n1,n2,n3                ! generic indices
   integer(IN)  :: kui,kvi                 ! field indices
   integer(IN)  :: kuo,kvo                 ! field indices
   integer(IN)  :: kin                     ! index index
   integer(IN)  :: nmin,nmax               ! indices of highest latitude in input
   integer(IN)  :: npts                    ! local number of points in an aV
   integer(IN)  :: num                     ! number of points at highest latitude
   integer(IN)  :: kloni                   ! longitude index on input domain
   integer(IN)  :: klati                   ! latitude index on input domain
   integer(IN)  :: klono                   ! longitude index on output domain
   integer(IN)  :: klato                   ! latitude index on output domain
   integer(IN)  :: index                   ! index value
   real(R8)     :: rindex                  ! index value
   real(R8)     :: latmax                  ! value of highest latitude
   real(R8)     :: olon,olat               ! output bundle lon/lat
   real(R8)     :: ilon,ilat               ! input bundle lon/lat
   real(R8)     :: npu,npv                 ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2           ! angles for trig functions
   real(R8),allocatable :: rarray(:)       ! temporary array
   real(R8),allocatable :: rarray2(:,:)    ! temporary array
   real(R8)     :: w1,w2,w3,w4             ! weights
   real(R8)     :: f1,f2,f3,f4             ! function values
   real(R8)     :: alpha,beta              ! used to generate weights
   real(R8)     :: rtmp                    ! real temporary
   real(R8),allocatable :: lData(:,:)      ! last lat local input bundle data
   integer      :: mpicom                  ! mpi communicator group 
   real(R8)   ,allocatable :: rfound(:)    ! temporary for copy
   integer(IN),allocatable :: ifound(:)    ! temporary for copy
   integer(IN)  :: cnt                     ! loop counter
   logical      :: found                   ! search for new interpolation
   integer(IN)  :: rcode                   ! error code
   integer(IN)  :: np1                     ! n+1 or tmp
   integer(IN)  :: mype

   !--- formats ---
   character(*),parameter :: subName = '(map_npFixNew4R) '
   character(*),parameter :: F00 = "('(map_npFixNew4R) ',8a)"
   character(*),parameter :: F01 = "('(map_npFixNew4R) ',a,i12)"
   !-------------------------------------------------------------------------------

   mpicom = mapper%mpicom
   call MPI_COMM_RANK(mpicom,mype,rcode)

   nmin = (mapper%ni)*(mapper%nj-1) + 1
   nmax = mapper%ni*mapper%nj
   num  = mapper%ni

   !---------------------------------------------------------------------------
   ! Initialization only
   !---------------------------------------------------------------------------

   if (trim(mode) == 'init') then

     if (loglevel > 0) write(logunit,F00) " compute bilinear weights & indicies for NP region."

!    tcx 3/19/08, don't use GlobGridNum, it's not set properly in models
!    kin   = mct_aVect_indexIA(domi%data,"GlobGridNum",perrWith=subName)
     klati = mct_aVect_indexRA(domi%data,"lat"        ,perrWith=subName)
     kloni = mct_aVect_indexRA(domi%data,"lon"        ,perrWith=subName)
     klato = mct_aVect_indexRA(domo%data,"lat"        ,perrWith=subName)
     klono = mct_aVect_indexRA(domo%data,"lon"        ,perrWith=subName)

     allocate(mapper%ilon1(num))
     allocate(mapper%ilon2(num))
     allocate(mapper%ilat1(num))
     allocate(mapper%ilat2(num))

     mapper%ilon1 = 0._r8
     mapper%ilon2 = 0._r8
     mapper%ilat1 = 0._r8
     mapper%ilat2 = 0._r8

     npts = mct_aVect_lSize(domi%data)
     call mct_gsMap_orderedPoints(gsMapi, mype, mapper%gindex)
     do m=1,npts
       if (mapper%gindex(m).ge.nmin) then               ! are on highest latitude
         n = mapper%gindex(m) - nmin + 1                ! n is 1->mapper%ni lon index on highest latitude
         rtmp = domi%data%rAttr(kloni,m)      ! rtmp is longitude value on highest latitude
         mapper%ilon1(n) = mod(rtmp+360._R8,360._R8) ! ilon1(n) is longitude val mapped from 0->360
         mapper%ilat1(n) = domi%data%rAttr(klati,m)  ! ilat1(n) values should all be the same (i.e. highest lat)
       endif
     enddo

     !--- all gather local data, MPI_SUM is low memory and simple
     !--- but is a performance penalty compared to gatherv and copy
     !--- or a fancy send/recv 

     allocate(rarray(num))
     rarray = mapper%ilat1
     call MPI_ALLREDUCE(rarray,mapper%ilat1,num,MPI_REAL8,MPI_SUM,mpicom,rcode)
     if (rcode.ne.0) then
       write(logunit,*) trim(subName),' ilat1 rcode error ',rcode
       call shr_sys_abort()
     endif
     rarray = mapper%ilon1
     call MPI_ALLREDUCE(rarray,mapper%ilon1,num,MPI_REAL8,MPI_SUM,mpicom,rcode)
     if (rcode.ne.0) then
       write(logunit,*) trim(subName),' ilon1 rcode error ',rcode
       call shr_sys_abort()
     endif

     do n = 1,num
       np1 = mod(n,num)+1
       mapper%ilat2(n) = mapper%ilat1(np1)
       mapper%ilon2(n) = mapper%ilon1(np1)
       if (mapper%ilon2(n) < mapper%ilon1(n)) mapper%ilon2(n) = mapper%ilon2(n) + 360._R8
     enddo

     do n = 1,num
       if (mapper%ilat1(n) /= mapper%ilat2(n)) then
          write(logunit,*) trim(subname),' ERROR: ilat1 ne ilat2 ',n,mapper%ilat1(n),mapper%ilat2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilat1 ne ilat2')
       endif
       if (mapper%ilon2(n) < mapper%ilon1(n)) then
          write(logunit,*) trim(subname),' ERROR: ilon2 lt ilon1 ',n,mapper%ilon1(n),mapper%ilon2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilon2 ilon1 error')
       endif
       ! tcraig check that lon diffs are reasonable 4x average dlon seems like reasonable limit
       if (mapper%ilon2(n) - mapper%ilon1(n) > (360.0_R8/(num*1.0_R8))*4.0) then
          write(logunit,*) trim(subname),' ERROR: ilon1,ilon2 ',n,mapper%ilon1(n),mapper%ilon2(n)
          call shr_sys_abort(trim(subname)//' ERROR: ilon2 ilon1 size diff')
       endif
     enddo

     latmax = maxval(mapper%ilat1)

     !--- compute weights and save them ---

     npts = mct_aVect_lSize(domo%data)
     allocate(mapper%mfound(npts),mapper%nfound(npts),mapper%alphafound(npts),mapper%betafound(npts))
     mapper%cntfound = 0

     do m = 1,npts
       olat = domo%data%rAttr(klato,m)
       if (olat >= latmax) then
         rtmp = domo%data%rAttr(klono,m)
         olon = mod(rtmp,360._R8)
         n = 1
         found = .false.
         do while (n <= num .and. .not.found )
           if (    olon         >= mapper%ilon1(n) .and. olon < mapper%ilon2(n) .or.   &
                   olon+360._R8 >= mapper%ilon1(n) .and. olon < mapper%ilon2(n)) then


!tcx ilat2==ilat1 so don't average
!--->        ilat = (mapper%ilat1(n) + mapper%ilat2(n)) * 0.5_R8
             ilat = mapper%ilat1(n)
             if (mapper%ilon2(n) == mapper%ilon1(n)) then
               alpha = 0.5_R8
             else if (    olon >= mapper%ilon1(n) .and. olon < mapper%ilon2(n)) then
               alpha = (olon - mapper%ilon1(n)) / (mapper%ilon2(n) - mapper%ilon1(n))
             else if (olon+360._R8>= mapper%ilon1(n) .and. olon < mapper%ilon2(n)) then
               alpha = (olon+360._R8 - mapper%ilon1(n)) / (mapper%ilon2(n) - mapper%ilon1(n))
             else
               write(logunit,*) subName,' ERROR: olon ',olon,mapper%ilon1(n),mapper%ilon2(n)
             endif
             if (ilat >= 90._R8) then
               beta  = 1.0_R8
             else
               beta  = (olat - ilat) / (90._R8 - ilat)
             endif

             mapper%cntfound = mapper%cntfound + 1
             mapper%mfound(mapper%cntfound) = m
             mapper%nfound(mapper%cntfound) = n
             mapper%alphafound(mapper%cntfound) = alpha
             mapper%betafound(mapper%cntfound) = beta
             found = .true.

           endif
           n = n + 1     ! normal increment
         enddo
         if ( .not.found ) then
           write(logunit,*) subName,' ERROR: found = false ',found,m,olon,olat
         endif
       endif
     end do

     allocate(ifound(npts))
     ifound(1:mapper%cntfound) = mapper%mfound(1:mapper%cntfound)
     deallocate(mapper%mfound)
     if (mapper%cntfound > 0) then
        allocate(mapper%mfound(mapper%cntfound))
        mapper%mfound(1:mapper%cntfound) = ifound(1:mapper%cntfound)
     endif

     ifound(1:mapper%cntfound) = mapper%nfound(1:mapper%cntfound)
     deallocate(mapper%nfound)
     if (mapper%cntfound > 0) then
        allocate(mapper%nfound(mapper%cntfound))
        mapper%nfound(1:mapper%cntfound) = ifound(1:mapper%cntfound)
     endif
     deallocate(ifound)

     allocate(rfound(npts))
     rfound(1:mapper%cntfound) = mapper%alphafound(1:mapper%cntfound)
     deallocate(mapper%alphafound)
     if (mapper%cntfound > 0) then
        allocate(mapper%alphafound(mapper%cntfound))
        mapper%alphafound(1:mapper%cntfound) = rfound(1:mapper%cntfound)
     endif

     rfound(1:mapper%cntfound) = mapper%betafound(1:mapper%cntfound)
     deallocate(mapper%betafound)
     if (mapper%cntfound > 0) then
        allocate(mapper%betafound(mapper%cntfound))
        mapper%betafound(1:mapper%cntfound) = rfound(1:mapper%cntfound)
     endif
     deallocate(rfound)

     call MPI_ALLREDUCE(mapper%cntfound,mapper%cntf_tot,1,MPI_INTEGER,MPI_SUM,mpicom,rcode)
     if (mype == 0) then
        write(logunit,F01) ' total npfix points found = ',mapper%cntf_tot
     endif

   elseif (trim(mode) == 'run') then

      !---------------------------------------------------------------------------
      ! Return if there is nothing to do; must be no points on any pes
      ! If there are any npfix points, all pes must continue to the np u,v calc
      !---------------------------------------------------------------------------

      if (mapper%cntf_tot < 1) then
         return
      endif

      !---------------------------------------------------------------------------
      ! Non-initialization, run-time fix
      !---------------------------------------------------------------------------

      !--- barrier not required but interesting for timing. ---
      !  call shr_mpi_barrier(mpicom,subName//" barrier")

      !--- extract index, u, v from buni ---
      kui   = mct_aVect_indexRA(buni,fld1,perrWith=subName)
      kvi   = mct_aVect_indexRA(buni,fld2,perrWith=subName)
      kuo   = mct_aVect_indexRA(buno,fld1,perrWith=subName)
      kvo   = mct_aVect_indexRA(buno,fld2,perrWith=subName)

      allocate(lData(3,num))
      lData = 0._R8
      npts = mct_aVect_lSize(buni)
      do n=1,npts
        if (mapper%gindex(n).ge.nmin) then
          m = mapper%gindex(n) - nmin + 1
          lData(1,m) = mapper%gindex(n)
          lData(2,m) = buni%rAttr(kui,n)
          lData(3,m) = buni%rAttr(kvi,n)
        endif
      enddo
   
      !--- all gather local data, MPI_SUM is low memory and simple
      !--- but is a performance penalty compared to gatherv and copy
      !--- KLUDGE - this should be looked at when it becomes a performance/memory
      !--- penalty   
   
      allocate(rarray2(3,num))
      rarray2=lData
      call MPI_ALLREDUCE(rarray2,lData,3*num,MPI_REAL8,MPI_SUM,mpicom,rcode)
      deallocate(rarray2)
   
      if (rcode.ne.0) then
        write(logunit,*) trim(subName),' rcode error ',rcode
        call shr_sys_abort()
      endif
   
      do n2=1,num
        if (lData(1,n2).lt.0.1_R8) then
          write(logunit,*) trim(subName),' error allreduce ',n2
        endif
      enddo
   
      !--- compute npu, npv (pole data) and initialize ilon,ilat arrays ---
   
      npu = 0._R8
      npv = 0._R8
      do n = 1,num
        theta1 = mapper%ilon1(n)*deg2rad
        npu = npu + cos(theta1)*lData(2,n) &
                  - sin(theta1)*lData(3,n)
        npv = npv + sin(theta1)*lData(2,n) &
                  + cos(theta1)*lData(3,n)
      enddo
      npu = npu / real(num,R8)
      npv = npv / real(num,R8)
   
      !--- compute updated pole vectors ---
   
   !DIR$ CONCURRENT
      do cnt = 1,mapper%cntfound
         m     = mapper%mfound(cnt)
         n     = mapper%nfound(cnt)
         np1   = mod(n,num)+1
         alpha = mapper%alphafound(cnt)
         beta  = mapper%betafound(cnt)
   
         w1 = (1.0_R8-alpha)*(1.0_R8-beta)
         w2 = (    alpha)*(1.0_R8-beta)
         w3 = (    alpha)*(    beta)
         w4 = (1.0_R8-alpha)*(    beta)
   
         theta1 = mapper%ilon1(n)*deg2rad
         theta2 = mapper%ilon2(n)*deg2rad
   
         f1 = lData(2,n)
         f2 = lData(2,np1)
         f3 =  cos(theta1)*npu + sin(theta1)*npv
         f4 =  cos(theta2)*npu + sin(theta2)*npv
         rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
         buno%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
   
         f1 = lData(3,n)
         f2 = lData(3,np1)
         f3 = -sin(theta1)*npu + cos(theta1)*npv
         f4 = -sin(theta2)*npu + cos(theta2)*npv
         rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
         buno%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
      enddo
   
      deallocate(lData)

   else

      call shr_sys_abort(trim(subname)//' ERROR: mode not valid '//trim(mode))

   endif

end subroutine map_npFixNew4R

!===============================================================================
 subroutine seq_map_gsmapcheck(gsmap1,gsmap2)

    implicit none
    type(mct_gsMap),intent(in) :: gsmap1
    type(mct_gsMap),intent(in) :: gsmap2

    integer(IN) :: s1, s2
    character(*),parameter :: subName = '(seq_map_gsmapcheck) '

    s1 = mct_gsMap_gsize(gsMap1)
    s2 = mct_gsMap_gsize(gsMap2)
    if (s1 /= s2) then
      write(logunit,*) trim(subname),'gsmap global sizes different ',s1,s2
      call shr_sys_abort(subName // "different gsmap size")
    endif

 end subroutine seq_map_gsmapcheck
!===============================================================================

end module seq_map_mod
