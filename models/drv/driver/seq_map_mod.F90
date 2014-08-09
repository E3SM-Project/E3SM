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
  use shr_mct_mod, only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use mct_mod
  use seq_comm_mct
#ifdef USE_ESMF_LIB
  use esmf
  use esmfshr_mod
  use seq_map_esmf
#endif
  use component_type_mod
  use seq_map_type_mod

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_map_init_rcfile     ! cpl pes
  public :: seq_map_init_rearrolap  ! cpl pes
  public :: seq_map_initvect        ! cpl pes
  public :: seq_map_map             ! cpl pes 
  public :: seq_map_mapvect         ! cpl pes
  public :: seq_map_readdata        ! cpl pes
#ifdef USE_ESMF_LIB
  public :: seq_map_register        ! cpl pes 
#endif

  interface seq_map_avNorm 
     module procedure seq_map_avNormArr
     module procedure seq_map_avNormAvF
  end interface

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

#ifdef USE_ESMF_LIB
  subroutine seq_map_register(petlist, comp, import_state, export_state)
  
    implicit none
  
    integer, pointer                  :: petlist(:)
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

 !===============================================================================

  subroutine seq_map_init_rcfile( mapper, comp_s, comp_d, &
       maprcfile, maprcname, maprctype, samegrid, string, esmf_map)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)        ,intent(inout),pointer :: mapper
    type(component_type) ,intent(inout)         :: comp_s 
    type(component_type) ,intent(inout)         :: comp_d 
    character(len=*)     ,intent(in)            :: maprcfile
    character(len=*)     ,intent(in)            :: maprcname
    character(len=*)     ,intent(in)            :: maprctype
    logical              ,intent(in)            :: samegrid
    character(len=*)     ,intent(in),optional   :: string
    logical              ,intent(in),optional   :: esmf_map
    !
    ! Local Variables
    !
    type(mct_gsmap), pointer    :: gsmap_s ! temporary pointers
    type(mct_gsmap), pointer    :: gsmap_d ! temporary pointers
    integer(IN)                 :: mpicom
    character(CX)               :: mapfile
    character(CL)               :: maptype
    integer(IN)                 :: mapid
    integer(IN)                 :: ssize,dsize
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
    character(len=*),parameter  :: subname = "(seq_map_init_rcfile) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    gsmap_s => component_get_gsmap_cx(comp_s)
    gsmap_d => component_get_gsmap_cx(comp_d)

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="copy")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          call seq_map_mapinit(mapper,mpicom)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
       endif

    elseif (samegrid) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          ! --- Initialize rearranger
          call seq_map_mapinit(mapper,mpicom)
          mapper%rearrange_only = .true.
          mapper%strategy = "rearrange"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
          call seq_map_gsmapcheck(gsmap_s, gsmap_d)
          call mct_rearr_init(gsmap_s, gsmap_d, mpicom, mapper%rearr)
       endif

    else

       ! --- Initialize Smatp
       call shr_mct_queryConfigFile(mpicom,maprcfile,maprcname,mapfile,maprctype,maptype)

       call seq_map_mapmatch(mapid,gsMap_s=gsMap_s,gsMap_d=gsMap_d,mapfile=mapfile,strategy=maptype)

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          call seq_map_mapinit(mapper,mpicom)
          mapper%mapfile = trim(mapfile)
          mapper%strategy= trim(maptype)
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)

          call shr_mct_sMatPInitnc(mapper%sMatp, mapper%gsMap_s, mapper%gsMap_d, trim(mapfile),trim(maptype),mpicom)
          if (present(esmf_map)) mapper%esmf_map = esmf_map

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

       endif  ! mapid >= 0
    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
          mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
    endif

  end subroutine seq_map_init_rcfile

  !=======================================================================

  subroutine seq_map_init_rearrolap(mapper, comp_s, comp_d, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)        ,intent(inout),pointer :: mapper
    type(component_type) ,intent(inout)         :: comp_s 
    type(component_type) ,intent(inout)         :: comp_d 
    character(len=*)     ,intent(in),optional   :: string
    !
    ! Local Variables
    !
    integer(IN)                :: mapid
    type(mct_gsmap), pointer   :: gsmap_s
    type(mct_gsmap), pointer   :: gsmap_d
    integer(IN)                :: mpicom
    character(len=*),parameter :: subname = "(seq_map_init_rearrolap) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
!       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    gsmap_s => component_get_gsmap_cx(comp_s)
    gsmap_d => component_get_gsmap_cx(comp_d)

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="copy")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          call seq_map_mapinit(mapper,mpicom)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
       endif

    else
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")

       if (mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          ! --- Initialize rearranger
          call seq_map_mapinit(mapper, mpicom)
          mapper%rearrange_only = .true.
          mapper%strategy = "rearrange"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
          call seq_map_gsmapcheck(gsmap_s, gsmap_d)
          call mct_rearr_init(gsmap_s, gsmap_d, mpicom, mapper%rearr)
       endif

    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
          mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
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
    type(seq_map)   ,intent(inout)       :: mapper
    type(mct_aVect) ,intent(in)          :: av_s
    type(mct_aVect) ,intent(inout)       :: av_d
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

  subroutine seq_map_initvect(mapper, type, comp_s, comp_d, string)

    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)        ,intent(inout)       :: mapper
    character(len=*)     ,intent(in)          :: type
    type(component_type) ,intent(inout)       :: comp_s 
    type(component_type) ,intent(inout)       :: comp_d 
    character(len=*)     ,intent(in),optional :: string
    !
    ! Local Variables
    !
    type(mct_gGrid), pointer   :: dom_s
    type(mct_gGrid), pointer   :: dom_d
    integer(IN)                :: klon, klat, lsize, n
    logical                    :: lnorm
    character(len=CL)          :: lstring
    character(len=*),parameter :: subname = "(seq_map_initvect) "
    !-----------------------------------------------------

    lstring = ' '
    if (present(string)) then
!       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
       lstring = trim(string)
    endif

    dom_s => component_get_dom_cx(comp_s)
    dom_d => component_get_dom_cx(comp_d)

    if (trim(type(1:6)) == 'cart3d') then
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

  subroutine seq_map_mapvect( mapper, type, av_s, av_d, fldu, fldv, norm, string )

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout)       :: mapper
    character(len=*),intent(in)          :: type
    type(mct_aVect) ,intent(in)          :: av_s
    type(mct_aVect) ,intent(inout)       :: av_d
    character(len=*),intent(in)          :: fldu
    character(len=*),intent(in)          :: fldv
    logical         ,intent(in),optional :: norm
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

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (trim(type(1:6)) == 'cart3d') then
       if (mapper%cart3d_init /= trim(seq_map_stron)) then
          call shr_sys_abort(trim(subname)//' ERROR: cart3d not initialized '//trim(lstring))
       endif
       call seq_map_cart3d(mapper, type, av_s, av_d, fldu, fldv, norm=lnorm, string=string)
    elseif (trim(type) == 'none') then
       call seq_map_map(mapper, av_s, av_d, fldlist=trim(fldu)//':'//trim(fldv), norm=lnorm)
    else
       write(logunit,*) subname,' ERROR: type unsupported ',trim(type)
       call shr_sys_abort(trim(subname)//' ERROR in type='//trim(type))
    end if

  end subroutine seq_map_mapvect

  !=======================================================================

  subroutine seq_map_cart3d( mapper, type, av_s, av_d, fldu, fldv, norm, string)

    implicit none
    !-----------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_map)   ,intent(inout)       :: mapper
    character(len=*),intent(in)          :: type
    type(mct_aVect) ,intent(in)          :: av_s
    type(mct_aVect) ,intent(inout)       :: av_d
    character(len=*),intent(in)          :: fldu
    character(len=*),intent(in)          :: fldv
    logical         ,intent(in),optional :: norm
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

       call seq_map_map(mapper, av3_s, av3_d, norm=lnorm)

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
    character(len=*),intent(in)             :: maprcfile
    character(len=*),intent(in)             :: maprcname
    integer(IN)     ,intent(in)             :: mpicom
    integer(IN)     ,intent(in)             :: ID
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
    type(seq_map)   , intent(inout)       :: mapper  ! mapper
    type(mct_aVect) , intent(in)          :: av_i    ! input 
    type(mct_aVect) , intent(inout)       :: av_o    ! output
    type(mct_aVect) , intent(in)          :: avf_i   ! extra src "weight"
    character(len=*), intent(in)          :: avfifld ! field name in avf_i
    character(len=*), intent(in),optional :: rList   ! fields list
    logical         , intent(in),optional :: norm    ! normalize at end
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

end module seq_map_mod
