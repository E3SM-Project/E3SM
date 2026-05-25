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
  use shr_kind_mod      ,only: CL => SHR_KIND_CL, CX => SHR_KIND_CX, CXX => SHR_KIND_CXX
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use shr_mpi_mod
  use mct_mod
  use seq_comm_mct
  use component_type_mod
  use seq_map_type_mod
  use seq_nlmap_mod
  use shr_moab_mod

  implicit none
  save
  private  ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: seq_map_init_rcfile     ! cpl pes
  public :: moab_map_init_rcfile    ! cpl pes
  public :: seq_map_initvect_moab   ! cpl pes
  public :: seq_map_map             ! cpl pes
  public :: seq_map_mapvect         ! cpl pes
  public :: seq_map_clean_moab      ! cpl pes




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

  subroutine seq_map_init_rcfile( mapper, comp_s, comp_d, &
       maprcfile, maprcname, maprctype, samegrid, string, esmf_map, no_match)

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
    logical              ,intent(in),optional   :: no_match
    !
    ! Local Variables
    !
    type(mct_gsmap), pointer    :: gsmap_s ! temporary pointers
    type(mct_gsmap), pointer    :: gsmap_d ! temporary pointers
    integer(IN)                 :: mpicom
    character(CX)               :: mapfile
    character(CL)               :: maptype
    integer(IN)                 :: mapid
    logical                     :: skip_match
    character(len=*),parameter  :: subname = "(seq_map_init_rcfile) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    skip_match = .false.
    if (present(no_match)) then
       if (no_match) skip_match = .true.
    endif
    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    gsmap_s => component_get_gsmap_cx(comp_s)
    gsmap_d => component_get_gsmap_cx(comp_d)

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then
       if(.not.skip_match) call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="copy")

       if (.not. skip_match .and. mapid > 0 ) then
          call seq_map_mappoint(mapid,mapper)
       else
          if(skip_match .and. seq_comm_iamroot(CPLID)) then
             write(logunit,'(A)') subname//' skip_match true, force new map'
          endif
          call seq_map_mapinit(mapper,mpicom)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
       endif

    elseif (samegrid) then
       if(.not.skip_match) call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")

       if (.not. skip_match .and. mapid > 0 ) then
          call seq_map_mappoint(mapid,mapper)
       else
          if(skip_match .and. seq_comm_iamroot(CPLID)) then
             write(logunit,'(A)') subname//'skip_match true, force new map'
          endif
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

       if(.not.skip_match) call seq_map_mapmatch(mapid,gsMap_s=gsMap_s,gsMap_d=gsMap_d,mapfile=mapfile,strategy=maptype)

       if (.not. skip_match .and. mapid > 0) then
          call seq_map_mappoint(mapid,mapper)
       else
          if(skip_match .and. seq_comm_iamroot(CPLID)) then
             write(logunit,'(A)') subname//'skip_match true, force new map'
          endif
          call seq_map_mapinit(mapper,mpicom)
          mapper%mapfile = trim(mapfile)
          mapper%strategy= trim(maptype)
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)

          call shr_mct_sMatPInitnc(mapper%sMatp, mapper%gsMap_s, mapper%gsMap_d, trim(mapfile),trim(maptype),mpicom)
          if (present(esmf_map)) mapper%esmf_map = esmf_map

          if (mapper%esmf_map) then
             call shr_sys_abort(subname//' ERROR: esmf SMM not supported')
          endif  ! esmf_map

       endif  ! mapid >= 0
    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
            mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
    endif

  end subroutine seq_map_init_rcfile


  subroutine moab_map_init_rcfile( mapper, discretization_type, &
                   maprcfile, maprcname, maprctype, samegrid, arearead, map_identifier, &
                   description_string, esmf_map, fallback_map_identifier )

    use iMOAB, only: iMOAB_LoadMapFile
    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_map)        ,intent(inout)         :: mapper  ! mapper being initialized (src_mbid, tgt_mbid, intx_mbid must be set)
    type(integer)        ,intent(in)            :: discretization_type ! 1 for SE, 2 for PC, 3 for FV; should be a member data
    ! type(component_type) ,intent(inout)         :: comp_s
    ! type(component_type) ,intent(inout)         :: comp_d
    character(len=*)     ,intent(in)            :: maprcfile
    character(len=*)     ,intent(in)            :: maprcname
    character(len=*)     ,intent(in)            :: maprctype
    logical              ,intent(in)            :: samegrid
    integer              ,intent(in)            :: arearead ! read or not area_a and area_b
    character(len=*)     ,intent(inout)         :: map_identifier !   /* "scalar", "flux", "custom" */
    character(len=*)     ,intent(in),optional   :: description_string
    logical              ,intent(in),optional   :: esmf_map
    character(len=*)     ,intent(in),optional   :: fallback_map_identifier
    !
    ! Local Variables
    !
    !type(mct_gsmap), pointer    :: gsmap_s ! temporary pointers
    !type(mct_gsmap), pointer    :: gsmap_d ! temporary pointers
    integer(IN)                 :: mpicom
    character(CX)               :: mapfile, nl_mapfile
    character(CX)               :: mapfile_term
    character(CL)               :: maptype
    integer(IN)                 :: mapid
    character(CL)               :: nlmap_id
    character(len=128)          :: nl_label
    logical                     :: nl_found
    integer                     :: ierr

    character(len=*),parameter  :: subname = "(moab_map_init_rcfile) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(description_string)) then
        write(logunit,'(A)') subname//' called for '//trim(description_string)
    endif

    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    ! --- Initialize Smatp
    call shr_mct_queryConfigFile(mpicom,maprcfile,maprcname,mapfile,maprctype,maptype)

    ! if this routine is called, there should be a mapfile present
    if (mapfile == 'idmap' .or. mapfile == 'idmap_ignore') then
      if (present(fallback_map_identifier)) then
         map_identifier = fallback_map_identifier
         write(logunit,*) subname,' got idmap. do not want to load backup identifier - ' // maprcname
         call shr_sys_abort(subname//' ERROR in not wanting to load backup identifier - ' // maprcname)
      else
         write(logunit,*) subname,' got idmap. expecting map file name for - ' // maprcname
         call shr_sys_abort(subname//' ERROR in finding map file - ' // maprcname)
      end if
    endif

    mapfile_term = trim(mapfile)//CHAR(0)
    if (seq_comm_iamroot(CPLID)) then
        write(logunit,*) subname,' reading map file with iMOAB: ', trim(mapfile_term)
    endif

    ierr = iMOAB_LoadMapFile( mapper%src_mbid, mapper%tgt_mbid, mapper%intx_mbid, discretization_type, &
                                 discretization_type, arearead, map_identifier, mapfile_term)
    if (ierr .ne. 0) then
       write(logunit,*) subname,' error in loading map file - ' // mapfile
       call shr_sys_abort(subname//' ERROR in loading map file - ' // mapfile)
    endif

    mapper%nl_available = .false.
    ! Look for an optional nonlinear (high-order) map paired with this one.
    nl_label = maprcname(1:len(maprcname)-1)//'_nonlinear:'
    call shr_mct_queryConfigFile(mpicom, maprcfile, trim(nl_label), nl_mapfile, &
          Label1Found=nl_found)
    if (nl_found) nl_found = nl_mapfile /= "idmap_ignore"

    if (nl_found) then
         mapper%nl_available = .true.
         nlmap_id = 'ho_'//map_identifier
         mapper%howeight_identifier = nlmap_id
         mapfile_term = trim(nl_mapfile)//CHAR(0)
         ierr = iMOAB_LoadMapFile( mapper%src_mbid, mapper%tgt_mbid, mapper%intx_mbid, discretization_type, &
                                 discretization_type, 0, nlmap_id, mapfile_term)
         if (ierr .ne. 0) then
           write(logunit,*) subname,' error in loading nlmap file - ' // nl_mapfile
           call shr_sys_abort(subname//' ERROR in loading nlmap file - ' // nl_mapfile)
         endif
    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
            mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       if (mapper%nl_available) then
          write(logunit,'(2A,I6,3A)') subname, &
               ' mapper counter, nlmap_id = ', &
               mapper%counter,' ',trim(nlmap_id)
       end if
       call shr_sys_flush(logunit)
    endif

  end subroutine moab_map_init_rcfile

  !=======================================================================

  !=======================================================================
  !
  ! seq_map_map - Maps attribute vector fields from source to destination grid.
  !
  ! This subroutine is the primary mapping interface for transferring field
  ! data between component grids within the coupler. It supports 2 mapping strategies:
  !
  ! 1. REARRANGE (mapper%rearrange_only or mapper%copy_only =.true.)
  !    - Used when grids have the SAME meshes but possibly DIFFERENT MPI distribution
  !    - Performs MPI communication to redistribute data
  !    - goes through MPI layer even if decomps are identical.
  !
  ! 2. MAP
  !    - Full interpolation/regridding between DIFFERENT grids
  !    - Uses sparse matrix multiplication with pre-computed or online weights
  !    - Supports optional normalization for conservation
  !
  ! MOAB uses iMOAB API calls to:
  ! - Transfer data between MOAB application instances
  ! - Apply projection weights for regridding
  ! - Handle normalization for conservative mapping
  !
  ! Arguments:
  !   mapper     - Mapping data structure containing weights and grid info
  !   av_s       - Source attribute vector (input fields)
  !   av_d       - Destination attribute vector (output fields)
  !   fldlist    - Optional: specific fields to map (default: all fields)
  !   norm       - Optional: enable normalization (default: .true.)
  !   avwts_s    - Optional: source weights for normalization
  !   avwtsfld_s - Optional: field name in avwts_s to use as weight
  !   string     - Optional: description string for logging
  !   msgtag     - Optional: MPI message tag for rearrangement
  !
  !=======================================================================

  subroutine seq_map_map( mapper, av_s, av_d, fldlist, norm, avwts_s, avwtsfld_s, &
       string, msgtag, omit_nonlinear  )

    use iso_c_binding
    use iMOAB, only: iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, iMOAB_SetDoubleTagStorage, &
      iMOAB_GetIntTagStorage, iMOAB_ApplyScalarProjectionWeights, &
      iMOAB_SendElementTag, iMOAB_ReceiveElementTag, iMOAB_FreeSenderBuffers
    use seq_comm_mct, only : num_moab_exports

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
    logical         ,intent(in),optional :: omit_nonlinear
    logical  :: valid_moab_context
    integer  :: ierr, nfields, lsize_src, lsize_tgt, arrsize_tgt, j, arrsize_src
    character(len=CXX) :: fldlist_moab
    character(len=CXX) :: fldlist_data    ! data fields only (no norm8wt) for nlmap CAAS path
    character(len=CXX) :: fldlist_caas    ! non-excluded data fields → dual-map CAAS
    character(len=CXX) :: fldlist_lo_only ! excluded data + norm8wt → low-order projection
    integer            :: ncaas_fields, nlo_fields  ! field counts for the two sub-lists
    character(len=CXX) :: tagname
    integer    :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
    type(mct_list) :: temp_list
    integer, dimension(:), allocatable  :: globalIds
    real(r8), dimension(:), allocatable  :: wghts
    real(kind=r8) , allocatable  :: targtags(:,:), targtags_ini(:,:)
    real(kind=r8)  :: factor
    integer  :: filter_type ! used for caas projection
    !
    ! Local Variables
    !
    logical :: lnorm  ! true if normalization is to be done
    logical :: mbnorm ! moab copy of lnorm
    logical :: use_nonlinear_map
    logical :: mbpresent ! moab logical for presence of norm weight string
    integer(IN),save :: ltag    ! message tag for rearrange
    character(len=*),parameter :: subname = "(seq_map_map) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    use_nonlinear_map = .false.
    if (mapper%nl_available) then
       use_nonlinear_map = .true.
       if (present(omit_nonlinear)) then
          if (omit_nonlinear) use_nonlinear_map = .false.
       end if
    end if

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    mbnorm = lnorm

    if (present(avwtsfld_s)) then
       mbpresent = .true.
    else
       mbpresent = .false.
    endif

    !mbnorm = .false.   ! uncomment both to turn off normalization for all maps
    !mbpresent = .false.

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

    !*** MOAB: Check if MOAB context is valid on this process
    !*** A valid MOAB context requires both source and target MOAB app IDs >= 0.
    !*** When valid, MOAB operations run IN PARALLEL with MCT operations,
    !*** maintaining data consistency between both representations.
    if ( mapper%src_mbid .lt. 0 .or. mapper%tgt_mbid .lt. 0 ) then
       valid_moab_context = .FALSE.
    else
       valid_moab_context = .TRUE.
    endif

    !*** MOAB: Build field list for MOAB operations
    !*** This section constructs a colon-delimited field list string that will
    !*** be passed to iMOAB functions. The list is built from either:
    !*** - The explicit fldlist argument, or
    !*** - All real attributes in the source attribute vector
    !*** If normalization is enabled (mbnorm=.true.), "norm8wt" is appended
    !*** to carry normalization weights through the mapping.
    if ( valid_moab_context ) then
       nfields = 1
       ! first get data from source tag and store in a temporary
       ! then set it back to target tag to mimic a copy
       if (present(fldlist)) then
          ! find the number of fields in the list
          ! Or should we decipher based on fldlist?
          call mct_list_init(temp_list, fldlist)
          nfields=mct_list_nitem (temp_list)
          call mct_list_clean(temp_list)
          fldlist_moab= trim(fldlist)
       else
          ! Extract character strings from attribute vector
          nfields = mct_aVect_nRAttr(av_s)
          fldlist_moab = ''
          if ( nfields /= 0 ) fldlist_moab = trim(mct_aVect_exportRList2c(av_s))
       endif

       ! Snapshot the data-only field list (no norm8wt). Used in the nlmap
       ! path so the dual-map CAAS only sees real data fields. norm8wt has to
       ! be mapped LOW-order separately (see Step 5 below) — putting it in the
       ! CAAS list mishandles it: the high-order map of constant 1.0 is not
       ! row-sum-1, so the CAAS-clipped result no longer matches MCT's
       ! mct_sMat_avMult-mapped target norm8wt that the post-divide expects.
       fldlist_data = trim(fldlist_moab)//C_NULL_CHAR

       if (mbnorm) then
          fldlist_moab = trim(fldlist_moab)//":norm8wt"//C_NULL_CHAR
          nfields=nfields + 1
       else
          fldlist_moab = trim(fldlist_moab)//C_NULL_CHAR
       endif


#ifdef MOABDEBUG
       if (seq_comm_iamroot(CPLID)) then
          write(logunit,*) subname, 'iMOAB mapper ',trim(mapper%mbname), ' iMOAB_mapper  nfields', &
                nfields,  ' fldlist_moab=', trim(fldlist_moab), ' moab step ', num_moab_exports
          call shr_sys_flush(logunit)
       endif
#endif
    endif ! valid_moab_context

    !=========================================================================
    ! MOAB PATH: Copy/Rearrange Operations
    ! For COPY and REARRANGE strategies, MOAB uses point-to-point
    ! communication between MOAB app instances instead of MCT rearranger.
    !=========================================================================
    if (mapper%copy_only .or. mapper%rearrange_only) then

       !*** MOAB: Send/Receive for Copy or Rearrange
       !*** For identical or same-grid mappings, MOAB transfers data directly
       !*** between source and target MOAB apps using iMOAB_SendElementTag
       !*** and iMOAB_ReceiveElementTag. This is analogous to MCT's copy or
       !*** rearrange but operates on MOAB mesh data structures.
       if ( valid_moab_context ) then
#ifdef MOABDEBUG
          if (seq_comm_iamroot(CPLID)) then
             write(logunit, *) subname,' iMOAB mapper rearrange or copy ', mapper%mbname, ' send/recv tags ', trim(fldlist_moab), &
               ' mbpresent=', mbpresent, ' mbnorm=', mbnorm, ' moab step:', num_moab_exports
             call shr_sys_flush(logunit)
          endif
#endif
          !*** MOAB: iMOAB_SendElementTag - Send field data from source MOAB app
          !*** Parameters:
          !***   src_mbid: Source MOAB application ID
          !***   fldlist_moab: Colon-delimited list of tag names to send
          !***   mpicom: MPI communicator
          !***   intx_context: Context ID for the intersection/target app
          ierr = iMOAB_SendElementTag( mapper%src_mbid, fldlist_moab, mapper%mpicom, mapper%intx_context );
          if (ierr .ne. 0) then
             write(logunit, *) subname,' iMOAB mapper ', mapper%mbname, ' error in sending tags ', trim(fldlist_moab), ierr
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//' ERROR in sending tags')
          endif
          !*** MOAB: iMOAB_ReceiveElementTag - Receive field data into target MOAB app
          !*** Parameters:
          !***   tgt_mbid: Target MOAB application ID
          !***   fldlist_moab: Colon-delimited list of tag names to receive
          !***   mpicom: MPI communicator
          !***   src_context: Context ID for the source app
          ierr = iMOAB_ReceiveElementTag( mapper%tgt_mbid, fldlist_moab, mapper%mpicom, mapper%src_context );
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in receiving tags iMOAB mapper ', mapper%mbname,  trim(fldlist_moab)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//' ERROR in receiving tags')
          endif
          !*** MOAB: iMOAB_FreeSenderBuffers - Release MPI send buffers
          !*** Must be called after receive completes to free memory.
          ierr = iMOAB_FreeSenderBuffers( mapper%src_mbid, mapper%intx_context )
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in freeing buffers ', trim(fldlist_moab)
             call shr_sys_abort(subname//' ERROR in freeing buffers') ! serious enough
          endif
       endif ! if (valid_moab_context)


    else
       !=========================================================================
       ! MOAB PATH: Full Mapping Operations (Strategy 3)
       ! For full interpolation/regridding between different grids.
       ! This is more complex than copy/rearrange and involves:
       ! 1. Pre-normalization: Multiply source fields by normalization weight
       ! 2. Send/Receive: Transfer data to intersection mesh
       ! 3. Apply weights: Apply projection/interpolation weights
       ! 4. Post-normalization: Divide by mapped normalization weight
       !=========================================================================

       !*** MOAB: Pre-normalization (source side)
       !*** For conservative mapping, we need to:
       !*** 1. Set norm8wt tag to 1.0 on all source cells
       !*** 2. If a weight field (avwtsfld_s) is specified, multiply all fields by it
       !*** 3. Map both the fields AND the weights
       !*** 4. Post-map: divide fields by mapped weights to restore proper scaling
       if ( valid_moab_context ) then
          !*** MOAB: Pre-normalization setup
          !*** When normalization is requested (mbnorm=.true.) or a weight field
          !*** is provided (mbpresent=.true.), we prepare source data for mapping.
          if (mbnorm .or. mbpresent) then
             !*** MOAB: iMOAB_GetMeshInfo - Get source mesh dimensions
             !*** nvise(1) returns the number of visible (owned) elements
             ierr  = iMOAB_GetMeshInfo ( mapper%src_mbid, nvert, nvise, nbl, nsurf, nvisBC );
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error getting mesh info for ', mapper%mbname
                call shr_sys_abort(subname//' ERROR getting mesh info') ! serious enough
             endif
             lsize_src = nvise(1) ! number of active cells

             !*** MOAB: Initialize norm8wt tag to 1.0 on all source cells
             !*** This tag will be mapped along with the data fields.
             !*** After mapping, dividing by the mapped norm8wt gives proper normalization.
             allocate(wghts(lsize_src))
             wghts = 1.0_r8
             tagname = "norm8wt"//C_NULL_CHAR
             !*** MOAB: iMOAB_SetDoubleTagStorage - Set tag values on mesh elements
             ierr = iMOAB_SetDoubleTagStorage (mapper%src_mbid, tagname, lsize_src , mapper%tag_entity_type, wghts)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error setting init value for mapping norm factor ',ierr,trim(tagname)
                call shr_sys_abort(subname//' ERROR setting norm init value') ! serious enough
             endif
#ifdef MOABDEBUG
             if (seq_comm_iamroot(CPLID)) then
                write(logunit, *) subname,' iMOAB mapper ', mapper%mbname, ' set norm8wt 1  on source with app id: ', &
                    mapper%src_mbid, ' moab step:', num_moab_exports
                call shr_sys_flush(logunit)
             endif
#endif
             !*** MOAB: Optional pre-multiplication by user-specified weights
             !*** If avwtsfld_s was provided, multiply all source fields by this weight.
             !*** This is used for flux-weighted mapping (e.g., area-weighted averages).
             if(mbpresent) then
                !*** MOAB: iMOAB_GetDoubleTagStorage - Get weight field values
                tagname = avwtsfld_s//C_NULL_CHAR
                ierr = iMOAB_GetDoubleTagStorage (mapper%src_mbid, tagname, lsize_src , mapper%tag_entity_type, wghts)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error getting value for mapping norm factor ', trim(tagname)
                   call shr_sys_abort(subname//' ERROR getting norm factor') ! serious enough
                endif

                !*** MOAB: Get all source field values and save original values
                !*** targtags_ini stores original values to restore after mapping
                allocate(targtags(lsize_src,nfields))
                allocate(targtags_ini(lsize_src,nfields))
                arrsize_src=lsize_src*(nfields)

                !*** MOAB: iMOAB_GetDoubleTagStorage - Get current field values
                ierr = iMOAB_GetDoubleTagStorage (mapper%src_mbid, fldlist_moab, arrsize_src , mapper%tag_entity_type, targtags)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error getting source tag values ', mapper%mbname, mapper%src_mbid, trim(fldlist_moab), arrsize_src, mapper%tag_entity_type
                   call shr_sys_abort(subname//' ERROR getting source tag values') ! serious enough
                endif

                targtags_ini = targtags
                !*** MOAB: Multiply all fields by the normalization weight
                !*** Since norm8wt starts at 1.0, after this multiplication:
                !*** - Data fields contain: field_value * weight
                !*** - norm8wt contains: 1.0 * weight = weight
                !*** After mapping and dividing by mapped norm8wt, we get proper normalization.
                do j = 1, lsize_src
                   targtags(j,:)= targtags(j,:)*wghts(j)
                enddo
#ifdef MOABDEBUG
                if (seq_comm_iamroot(CPLID)) then
                   write(logunit, *) subname,' iMOAB projection mapper: ', mapper%mbname, ' normalize nfields=', &
                      nfields, ' arrsize_src on root:', arrsize_src, ' shape(targtags_ini)=', shape(targtags_ini), &
                      ' moab step:', num_moab_exports
                   call shr_sys_flush(logunit)
                endif
#endif
                !*** MOAB: iMOAB_SetDoubleTagStorage - Store weighted values back to source
                !*** These weighted values will be sent and mapped to the target grid.
                ierr = iMOAB_SetDoubleTagStorage (mapper%src_mbid, fldlist_moab, arrsize_src , mapper%tag_entity_type, targtags)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error setting normed source tag values ', mapper%mbname
                   call shr_sys_abort(subname//' ERROR setting normed source tag values') ! serious enough
                endif
                deallocate(targtags)
             endif ! end multiplication by norm factor
             deallocate(wghts)
          endif  ! end NORMALIZATION

          !*** MOAB: Send source data to intersection/coverage mesh
          !*** For full mapping, data goes to the intersection app (intx_context)
          !*** where projection weights will be applied.
          ierr = iMOAB_SendElementTag( mapper%src_mbid, fldlist_moab, mapper%mpicom, mapper%intx_context )
          if (ierr .ne. 0) then
             write(logunit, *) subname,' iMOAB mapper error in sending tags ', mapper%mbname,  trim(fldlist_moab)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//' ERROR in sending tags')
          endif
       endif

       !*** MOAB: Receive data into intersection mesh for weight application
       !*** The data is received by the intersection app (intx_mbid), which holds
       !*** the projection weights. Note: for true intersection cases, tgt_mbid
       !*** may be the same as intx_mbid.
       if ( valid_moab_context ) then
#ifdef MOABDEBUG
          if (seq_comm_iamroot(CPLID)) then
             write(logunit, *) subname,' iMOAB mapper receiving tags with intx and intx_mbid: ', &
                mapper%mbname, trim(fldlist_moab), ' moab step:', num_moab_exports
          endif
#endif
          ierr = iMOAB_ReceiveElementTag( mapper%intx_mbid, fldlist_moab, mapper%mpicom, mapper%src_context )
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in receiving tags ', mapper%mbname, 'recv:',  mapper%intx_mbid, trim(fldlist_moab)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//' ERROR in receiving tags')
             !valid_moab_context = .false. ! do not attempt to project
          endif
          !*** MOAB: iMOAB_FreeSenderBuffers - Release MPI send buffers
          ierr = iMOAB_FreeSenderBuffers( mapper%src_mbid, mapper%intx_context )
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in freeing buffers ', trim(fldlist_moab)
             call shr_sys_abort(subname//' ERROR in freeing buffers') ! serious enough
          endif
       endif

       !*** MOAB: Apply projection weights to interpolate source data to target grid
       !*** This is the core remapping operation using iMOAB_ApplyScalarProjectionWeights.
       !*** The weights were loaded earlier via moab_map_init_rcfile.
       if ( valid_moab_context ) then

#ifdef MOABDEBUG
          if (seq_comm_iamroot(CPLID)) then
             write(logunit, *) subname,' iMOAB projection mapper: ',trim(mapper%mbname), ' between ', mapper%src_mbid, ' and ',  mapper%tgt_mbid, trim(fldlist_moab), &
                ' moab step:', num_moab_exports
             call shr_sys_flush(logunit)
          endif
#endif
          !*** MOAB: iMOAB_ApplyScalarProjectionWeights - Apply interpolation/projection
          !*** Parameters:
          !***   intx_mbid: Intersection mesh app ID (holds the weights)
          !***   filter_type: 0=no filter, other values for CAAS projection
          !***   weight_identifier: Name of the weight matrix (e.g., "scalar", "flux")
          !***   fldlist_moab: Input and output field names (can be different)
          if(.not.use_nonlinear_map) then
             filter_type = 0 ! CAAS_NONE: plain low-order projection
             ierr = iMOAB_ApplyScalarProjectionWeights ( mapper%intx_mbid, filter_type, mapper%weight_identifier, &
               fldlist_moab, fldlist_moab)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error in applying weights '
                call shr_sys_abort(subname//' ERROR in applying weights')
             endif
          else
             ! Dual-map nonlinear remapping path. To match MCT's
             ! seq_nlmap_avNormArr exactly we split the data fields by
             ! whether they appear in the namelist nlmaps_exclude_fields list:
             !
             !   (a) NON-EXCLUDED data fields (`fldlist_caas`): high-order
             !       interpolation with low-order CAAS bounds via the dual-map
             !       call. filter_type must be != CAAS_NONE for
             !       ApplyWeightsWithDualMap to actually run (otherwise iMOAB
             !       falls through to plain high-order without bounds —
             !       produces OOB values that crash icepack).
             !   (b) EXCLUDED data fields (`fldlist_lo_only`): low-order map
             !       only, matching MCT seq_nlmap_avNormArr's behavior at
             !       lines 691-698 of driver-mct/main/seq_nlmap_mod.F90 — for
             !       these fields avp_o keeps the low-order mct_sMat_avMult
             !       result and is NOT overwritten with the nonlinear-fixer
             !       output. Examples: Faxa_rainc/rainl/snowc/snowl.
             !   (c) `norm8wt` (when mbnorm): low-order map only, joined into
             !       the same low-order pass as (b) since both use the same
             !       weight matrix. Gives target norm8wt = low-order row sum.
             !       The post-norm divide below then divides data by this
             !       low-order row-sum, restoring correct normalization at
             !       partial-coverage cells (matches MCT outer
             !       seq_map_avNormArr post-divide).
             call build_nlmap_sublists( fldlist_data, mbnorm, &
                                        fldlist_caas, ncaas_fields, &
                                        fldlist_lo_only, nlo_fields )

             if (ncaas_fields > 0) then
                filter_type = 2 ! CAAS_LOCAL
                ierr = iMOAB_ApplyScalarProjectionWeights ( mapper%intx_mbid, filter_type, &
                       mapper%howeight_identifier, fldlist_caas, fldlist_caas )
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error in applying weights (data fields, dual-map CAAS)'
                   call shr_sys_abort(subname//' ERROR in applying weights (data fields, dual-map CAAS)')
                endif
             end if

             if (nlo_fields > 0) then
                ! Excluded data fields + (optionally) norm8wt — plain low-order.
                ! iMOAB takes the plain ApplyWeights branch (no dual-map CAAS)
                filter_type = 0
                ierr = iMOAB_ApplyScalarProjectionWeights ( mapper%intx_mbid, filter_type, &
                       mapper%weight_identifier, fldlist_lo_only, fldlist_lo_only )
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error in applying weights (excluded fields + norm8wt, low-order)'
                   call shr_sys_abort(subname//' ERROR in applying weights (excluded fields + norm8wt, low-order)')
                endif
             end if
          endif

          !*** MOAB: Post-normalization (target side)
          !*** After mapping, we need to divide by the mapped normalization weight.
          !*** This completes the conservative mapping: mapped_field / mapped_weight
          !*** gives the properly normalized field value on the target grid.
          if (mbnorm) then
             !*** MOAB: Get target mesh dimensions
             ierr  = iMOAB_GetMeshInfo ( mapper%tgt_mbid, nvert, nvise, nbl, nsurf, nvisBC );
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error getting mesh info for target ', mapper%mbname
                call shr_sys_abort(subname//' ERROR getting mesh info') ! serious enough
             endif

             lsize_tgt = nvise(1) ! number of active cells
             tagname = "norm8wt"//C_NULL_CHAR
             allocate(wghts(lsize_tgt))

             !*** MOAB: Get mapped normalization weights on target grid
             ierr = iMOAB_GetDoubleTagStorage (mapper%tgt_mbid, tagname, lsize_tgt , mapper%tag_entity_type, wghts)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error getting value for mapping norm factor post-map ', ierr, trim(tagname)
                call shr_sys_abort(subname//' ERROR getting norm factor') ! serious enough
             endif

             !*** MOAB: Get mapped field values on target grid
             allocate(targtags(lsize_tgt,nfields))
             arrsize_tgt=lsize_tgt*(nfields)
             ierr = iMOAB_GetDoubleTagStorage (mapper%tgt_mbid, fldlist_moab, arrsize_tgt , mapper%tag_entity_type, targtags)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error getting destination tag values ', mapper%mbname
                call shr_sys_abort(subname//' ERROR getting source tag values') ! serious enough
             endif

             !*** MOAB: Divide mapped fields by mapped normalization weight
             !*** For each target cell:
             !***   factor = 1/norm8wt if norm8wt != 0, else factor = norm8wt
             !***   final_value = mapped_value * factor
             !*** This gives the properly normalized field values.
             ! TODO:  add some check for wghts < puny
             do j = 1, lsize_tgt
                factor = wghts(j)
                if (wghts(j) .ne. 0) factor = 1.0_r8/wghts(j) ! should we compare to a small value instead ?
                targtags(j,:)= targtags(j,:)*factor
             enddo

             !*** MOAB: Store normalized field values back to target mesh
             ierr = iMOAB_SetDoubleTagStorage (mapper%tgt_mbid, fldlist_moab, arrsize_tgt , mapper%tag_entity_type, targtags)
             if (ierr .ne. 0) then
                write(logunit,*) subname,' error getting destination tag values ', mapper%mbname
                call shr_sys_abort(subname//' ERROR getting source tag values') ! serious enough
             endif

             deallocate(wghts, targtags)
             !*** MOAB: Restore original source field values
             !*** If we multiplied source fields by weights before mapping,
             !*** restore the original values so the source mesh is unchanged.
             if (mbpresent) then
#ifdef MOABDEBUG
                if (seq_comm_iamroot(CPLID)) then
                   write(logunit, *) subname,' iMOAB projection mapper: ', mapper%mbname, ' shape(targtags_ini)=', shape(targtags_ini)
                   call shr_sys_flush(logunit)
                endif
#endif
                !*** MOAB: Restore original values to source mesh
                ierr = iMOAB_SetDoubleTagStorage (mapper%src_mbid, fldlist_moab, arrsize_src , mapper%tag_entity_type, targtags_ini)
                if (ierr .ne. 0) then
                   write(logunit,*) subname,' error setting source tag values ', mapper%mbname
                   call shr_sys_abort(subname//' ERROR setting source tag values') ! serious enough
                endif
                deallocate(targtags_ini)
             endif

          endif ! end normalization

       endif

    endif ! end of mapping type if else

  end subroutine seq_map_map

  !=======================================================================

  subroutine seq_map_initvect_moab(mapper, type, src_mbid, tgt_mbid, string)
    use iMOAB, only: iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, iMOAB_DefineTagStorage
    use iso_c_binding, only: c_char, C_NULL_CHAR

    !-----------------------------------------------------
    !
    ! Purpose: Initialize vector mapping for MOAB data types
    ! This version uses iMOAB_GetDoubleTagStorage to retrieve
    ! lat/lon coordinates from MOAB mesh handles and computes
    ! trigonometric functions needed for spherical vector transformations.
    !
    ! The computed values are stored in MOAB-specific arrays:
    ! - slon_s_moab, clon_s_moab, slat_s_moab, clat_s_moab for source
    ! - slon_d_moab, clon_d_moab, slat_d_moab, clat_d_moab for target
    !
    ! These arrays are used by cart3d for vector field mapping
    ! between spherical coordinate systems.
    !
    ! Usage:
    !   call seq_map_initvect_moab(mapper, 'cart3d', src_app_id, tgt_app_id, 'my_mapper')
    !
    ! Arguments
    !
    type(seq_map)        ,intent(inout)       :: mapper
    character(len=*)     ,intent(in)          :: type     ! mapping type ('cart3d' for vector mapping)
    integer              ,intent(in)          :: src_mbid  ! source MOAB app ID
    integer              ,intent(in)          :: tgt_mbid  ! target MOAB app ID
    character(len=*)     ,intent(in),optional :: string   ! optional description string
    !
    ! Local Variables
    !
    integer(IN)                :: lsize_s, lsize_d, n, ierr,tagindex
    integer(IN)                :: num_verts, num_elems, num_edges, dimension
    character(len=CL)          :: lstring
    character(len=*),parameter :: subname = "(seq_map_initvect_moab) "
    character(64)              :: tagname
    integer                    :: ent_type
    real(R8), allocatable      :: lon_data(:), lat_data(:)
    !-----------------------------------------------------

    lstring = ' '
    if (present(string)) then
       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
       lstring = trim(string)
    endif

    if (trim(type(1:6)) == 'cart3d') then
       if (mapper%cart3d_init == trim(seq_map_stron)) return

       lsize_s = mbGetnCells(src_mbid)

       if (associated(mapper%slon_s_moab)) deallocate(mapper%slon_s_moab)
       if (associated(mapper%clon_s_moab)) deallocate(mapper%clon_s_moab)
       if (associated(mapper%slat_s_moab)) deallocate(mapper%slat_s_moab)
       if (associated(mapper%clat_s_moab)) deallocate(mapper%clat_s_moab)

       allocate(mapper%slon_s_moab(lsize_s), mapper%clon_s_moab(lsize_s), &
                mapper%slat_s_moab(lsize_s), mapper%clat_s_moab(lsize_s))
       allocate(lon_data(lsize_s), lat_data(lsize_s))

       ! Get longitude data from MOAB
       tagname = 'lon'//C_NULL_CHAR
       call mbGetCellTagVals(src_mbid, tagname,lon_data,lsize_s)

       ! Get latitude data from MOAB
       tagname = 'lat'//C_NULL_CHAR
       call mbGetCellTagVals(src_mbid, tagname,lat_data,lsize_s)

       ! Compute trigonometric functions for source
       do n = 1, lsize_s
          mapper%slon_s_moab(n) = sin(lon_data(n) * deg2rad)
          mapper%clon_s_moab(n) = cos(lon_data(n) * deg2rad)
          mapper%slat_s_moab(n) = sin(lat_data(n) * deg2rad)
          mapper%clat_s_moab(n) = cos(lat_data(n) * deg2rad)
       enddo

       deallocate(lon_data, lat_data)

       ! Get target mesh info and coordinates
       lsize_d = mbGetnCells(tgt_mbid)

       if (associated(mapper%slon_d_moab)) deallocate(mapper%slon_d_moab)
       if (associated(mapper%clon_d_moab)) deallocate(mapper%clon_d_moab)
       if (associated(mapper%slat_d_moab)) deallocate(mapper%slat_d_moab)
       if (associated(mapper%clat_d_moab)) deallocate(mapper%clat_d_moab)

       allocate(mapper%slon_d_moab(lsize_d), mapper%clon_d_moab(lsize_d), &
                mapper%slat_d_moab(lsize_d), mapper%clat_d_moab(lsize_d))
       allocate(lon_data(lsize_d), lat_data(lsize_d))

       ! Get longitude data from MOAB
       tagname = 'lon'//C_NULL_CHAR
       call mbGetCellTagVals(tgt_mbid, tagname,lon_data,lsize_d)

       ! Get latitude data from MOAB
       tagname = 'lat'//C_NULL_CHAR
       call mbGetCellTagVals(tgt_mbid, tagname,lat_data,lsize_d)

       ! Compute trigonometric functions for target
       do n = 1, lsize_d
          mapper%slon_d_moab(n) = sin(lon_data(n) * deg2rad)
          mapper%clon_d_moab(n) = cos(lon_data(n) * deg2rad)
          mapper%slat_d_moab(n) = sin(lat_data(n) * deg2rad)
          mapper%clat_d_moab(n) = cos(lat_data(n) * deg2rad)
       enddo

       deallocate(lon_data, lat_data)

       ! Define ux, uy, uz tags in source and destination MOAB apps
       tagname = 'ux:uy:uz'//C_NULL_CHAR
       ierr = iMOAB_DefineTagStorage(src_mbid, tagname, 1, 1, tagindex)
       if (ierr /= 0) then
          write(logunit,*) subname,' ERROR: Failed to define ux:uy:uz tags in source'
          call shr_sys_abort(trim(subname)//' ERROR defining source tags')
       endif

       ierr = iMOAB_DefineTagStorage(tgt_mbid, tagname, 1, 1, tagindex)
       if (ierr /= 0) then
          write(logunit,*) subname,' ERROR: Failed to define ux:uy:uz tags in target'
          call shr_sys_abort(trim(subname)//' ERROR defining target tags')
       endif
       mapper%cart3d_init = trim(seq_map_stron)
    endif

  end subroutine seq_map_initvect_moab


  subroutine seq_map_clean_moab(mapper)

    !-----------------------------------------------------
    !
    ! Purpose: Clean up MOAB-specific coordinate arrays
    !
    ! Deallocates and nullifies all MOAB coordinate arrays (slon_*_moab, etc.)
    ! and resets the cart3d_init flag. Should be called when the mapper
    ! is no longer needed to prevent memory leaks.
    !
    ! Usage:
    !   call seq_map_clean_moab(mapper)
    !
    ! Arguments
    !
    type(seq_map), intent(inout) :: mapper  ! mapper to clean up
    !
    ! Local Variables
    !
    character(len=*),parameter :: subname = "(seq_map_clean_moab) "
    !-----------------------------------------------------

    ! Deallocate MOAB coordinate arrays if allocated
    if (associated(mapper%slon_s_moab)) then
       deallocate(mapper%slon_s_moab)
       nullify(mapper%slon_s_moab)
    endif
    if (associated(mapper%clon_s_moab)) then
       deallocate(mapper%clon_s_moab)
       nullify(mapper%clon_s_moab)
    endif
    if (associated(mapper%slat_s_moab)) then
       deallocate(mapper%slat_s_moab)
       nullify(mapper%slat_s_moab)
    endif
    if (associated(mapper%clat_s_moab)) then
       deallocate(mapper%clat_s_moab)
       nullify(mapper%clat_s_moab)
    endif
    if (associated(mapper%slon_d_moab)) then
       deallocate(mapper%slon_d_moab)
       nullify(mapper%slon_d_moab)
    endif
    if (associated(mapper%clon_d_moab)) then
       deallocate(mapper%clon_d_moab)
       nullify(mapper%clon_d_moab)
    endif
    if (associated(mapper%slat_d_moab)) then
       deallocate(mapper%slat_d_moab)
       nullify(mapper%slat_d_moab)
    endif
    if (associated(mapper%clat_d_moab)) then
       deallocate(mapper%clat_d_moab)
       nullify(mapper%clat_d_moab)
    endif

    mapper%cart3d_init = trim(seq_map_stroff)

  end subroutine seq_map_clean_moab

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
       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
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

    use shr_moab_mod, only: mbGetCellTagVals, mbSetCellTagVals, mbGetnCells
    use iMOAB, only: iMOAB_DefineTagStorage
    use iso_c_binding
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
    integer           :: lsize,tagindex
    logical           :: lnorm
    integer           :: ku,kv,kux,kuy,kuz,n
    real(r8)          :: ue,un,ur,ux,uy,uz,speed
    real(r8)          :: urmaxl,urmax,uravgl,uravg,spavgl,spavg
    type(mct_aVect)   :: av3_s, av3_d
    integer(in)       :: mpicom,my_task,ierr,urcnt,urcntl
    ! MOAB-specific variables
    integer           :: lsize_s_moab, lsize_d_moab
    real(r8), allocatable :: u_vals_moab(:), v_vals_moab(:)
    real(r8), allocatable :: cart_vals_moab(:,:)  ! 2D array: (ux,uy,uz) x npoints
    logical           :: use_moab_data
    character(len=*),parameter :: subname = "(seq_map_cart3d) "
    character(64)              :: tagname

    lnorm = .true.
    if (present(norm)) then
       lnorm=norm
    endif

    mpicom = mapper%mpicom

    ! Check if we should use MOAB data (when src_mbid and tgt_mbid are set)
    use_moab_data = (mapper%src_mbid >= 0 .and. mapper%tgt_mbid >= 0)

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
    endif

    if (use_moab_data) then
       ! MOAB data path: get source values from MOAB
       lsize_s_moab = mbGetnCells(mapper%src_mbid)
       allocate(u_vals_moab(lsize_s_moab), v_vals_moab(lsize_s_moab))
       allocate(cart_vals_moab(lsize_s_moab,3))

       ! Get source u and v values from MOAB
       call mbGetCellTagVals(mapper%src_mbid, trim(fldu), u_vals_moab, lsize_s_moab)
       call mbGetCellTagVals(mapper%src_mbid, trim(fldv), v_vals_moab, lsize_s_moab)

       ! Convert source spherical to cartesian using MOAB coordinate arrays
       do n = 1, lsize_s_moab
          ur = 0.0_r8
          ue = u_vals_moab(n)
          un = v_vals_moab(n)
          ux = mapper%clon_s_moab(n)*mapper%clat_s_moab(n)*ur - &
            mapper%clon_s_moab(n)*mapper%slat_s_moab(n)*un - &
               mapper%slon_s_moab(n)*ue
          uy = mapper%slon_s_moab(n)*mapper%clat_s_moab(n)*ur - &
               mapper%slon_s_moab(n)*mapper%slat_s_moab(n)*un + &
               mapper%clon_s_moab(n)*ue
          uz = mapper%slat_s_moab(n)*ur + &
               mapper%clat_s_moab(n)*un
          cart_vals_moab(n,1) = ux
          cart_vals_moab(n,2) = uy
          cart_vals_moab(n,3) = uz
       enddo
       call mbSetCellTagVals(mapper%src_mbid, "ux:uy:uz"//C_NULL_CHAR, cart_vals_moab, lsize_s_moab*3)

    endif

    call seq_map_map(mapper, av3_s, av3_d, norm=lnorm)

    if (ku /= 0 .and. kv /= 0) then
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
       call mct_avect_clean(av3_s)
       call mct_avect_clean(av3_d)
    endif

    if (use_moab_data) then
       ! MOAB data path: convert cartesian back to spherical and set in MOAB
       lsize_d_moab = mbGetnCells(mapper%tgt_mbid)
       deallocate(u_vals_moab, v_vals_moab, cart_vals_moab)
       allocate(u_vals_moab(lsize_d_moab), v_vals_moab(lsize_d_moab))
       allocate(cart_vals_moab(lsize_d_moab,3))

       ! Get mapped values
       call mbGetCellTagVals(mapper%tgt_mbid, "ux:uy:uz"//C_NULL_CHAR, cart_vals_moab, lsize_d_moab*3)


       ! Convert back from cartesian to spherical on target grid using MOAB coordinates
       do n = 1, lsize_d_moab
          ux = cart_vals_moab(n,1)
          uy = cart_vals_moab(n,2)
          uz = cart_vals_moab(n,3)
          ue = -mapper%slon_d_moab(n)          *ux + &
               mapper%clon_d_moab(n)          *uy
          un = -mapper%clon_d_moab(n)*mapper%slat_d_moab(n)*ux - &
               mapper%slon_d_moab(n)*mapper%slat_d_moab(n)*uy + &
               mapper%clat_d_moab(n)*uz
          ur =  mapper%clon_d_moab(n)*mapper%clat_d_moab(n)*ux + &
               mapper%slon_d_moab(n)*mapper%clat_d_moab(n)*uy - &
               mapper%slat_d_moab(n)*uz
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
                u_vals_moab(n) = sign(speed,ue)
                v_vals_moab(n) = 0.0_r8
             else if (ue == 0.0_R8) then
                u_vals_moab(n) = 0.0_r8
                v_vals_moab(n) = sign(speed,un)
             else
                u_vals_moab(n) = sign(speed/sqrt(1.0_r8 + ((un*un)/(ue*ue))),ue)
                v_vals_moab(n) = sign(speed/sqrt(1.0_r8 + ((ue*ue)/(un*un))),un)
             endif
          else
             !--- this ignores ur ---
             u_vals_moab(n) = ue
             v_vals_moab(n) = un
          endif
       enddo

       ! Set final u,v values in destination MOAB app
       call mbSetCellTagVals(mapper%tgt_mbid, trim(fldu), u_vals_moab, lsize_d_moab)
       call mbSetCellTagVals(mapper%tgt_mbid, trim(fldv), v_vals_moab, lsize_d_moab)

       ! Clean up MOAB arrays
       deallocate(u_vals_moab, v_vals_moab, cart_vals_moab)
    endif


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

  end subroutine seq_map_cart3d

  !=======================================================================
  ! build_nlmap_sublists -- split a colon-separated tag list into two
  ! sub-lists for the dual-map nlmap path:
  !
  !   fldlist_caas    : data fields NOT in the nlmaps_exclude_fields list,
  !                     which go through the high-order CAAS dual-map call.
  !   fldlist_lo_only : data fields IN the exclude list, plus the special
  !                     "norm8wt" tag when mbnorm is true. These are mapped
  !                     with the LOW-order weight matrix only — same as MCT
  !                     seq_nlmap_avNormArr's behavior of leaving avp_o at
  !                     the mct_sMat_avMult result for excluded fields, and
  !                     of mapping the norm8wt column via the low-order map.
  !
  ! Both output strings are colon-separated, terminated with C_NULL_CHAR
  ! (so they can be passed straight to iMOAB). The corresponding field
  ! counts are also returned. The input fldlist must be C_NULL-terminated.
  !=======================================================================
  subroutine build_nlmap_sublists(fldlist_in, mbnorm, &
                                  fldlist_caas, ncaas, &
                                  fldlist_lo_only, nlo)
    use iso_c_binding, only : C_NULL_CHAR
    use seq_nlmap_mod, only : seq_nlmap_field_is_excluded

    character(len=*), intent(in)  :: fldlist_in
    logical,          intent(in)  :: mbnorm
    character(len=*), intent(out) :: fldlist_caas
    integer,          intent(out) :: ncaas
    character(len=*), intent(out) :: fldlist_lo_only
    integer,          intent(out) :: nlo

    integer :: i, n, start, ipos
    character(len=128) :: name

    fldlist_caas    = ''
    fldlist_lo_only = ''
    ncaas = 0
    nlo   = 0

    ! Trim trailing C_NULL_CHAR(s) before parsing.
    n = len_trim(fldlist_in)
    do while (n > 0)
       if (fldlist_in(n:n) /= C_NULL_CHAR) exit
       n = n - 1
    end do
    if (n <= 0) then
       fldlist_caas    = C_NULL_CHAR
       fldlist_lo_only = C_NULL_CHAR
       if (mbnorm) then
          fldlist_lo_only = 'norm8wt'//C_NULL_CHAR
          nlo = 1
       end if
       return
    end if

    ! Defense-in-depth: upper-bound buffer-size check. If every input field
    ! landed in one output list, that list would hold the entire input string
    ! (n chars including the colons between fields) plus the final
    ! C_NULL_CHAR -> n+1 chars. fldlist_lo_only additionally may need
    ! ':norm8wt' (8 chars) appended when mbnorm is true. If either output
    ! buffer is smaller than that worst case, abort with a clear message
    ! instead of silently truncating (which would manifest as an opaque
    ! "tag not found" error downstream in iMOAB).
    if (len(fldlist_caas) < n + 1) then
       call shr_sys_abort('(build_nlmap_sublists) ERROR: fldlist_caas buffer too small for input field list')
    end if
    if (len(fldlist_lo_only) < n + 1 + merge(8, 0, mbnorm)) then
       call shr_sys_abort('(build_nlmap_sublists) ERROR: fldlist_lo_only buffer too small for input field list')
    end if

    ! Walk the colon-separated list.
    start = 1
    do i = 1, n+1
       if (i == n+1 .or. fldlist_in(i:i) == ':') then
          if (i > start) then
             name = fldlist_in(start:i-1)
             if (seq_nlmap_field_is_excluded(name)) then
                if (nlo > 0) then
                   ipos = len_trim(fldlist_lo_only)
                   fldlist_lo_only(ipos+1:ipos+1) = ':'
                   fldlist_lo_only(ipos+2:) = trim(name)
                else
                   fldlist_lo_only = trim(name)
                end if
                nlo = nlo + 1
             else
                if (ncaas > 0) then
                   ipos = len_trim(fldlist_caas)
                   fldlist_caas(ipos+1:ipos+1) = ':'
                   fldlist_caas(ipos+2:) = trim(name)
                else
                   fldlist_caas = trim(name)
                end if
                ncaas = ncaas + 1
             end if
          end if
          start = i + 1
       end if
    end do

    ! Append norm8wt to the low-order list when normalization is requested.
    if (mbnorm) then
       if (nlo > 0) then
          ipos = len_trim(fldlist_lo_only)
          fldlist_lo_only(ipos+1:ipos+1) = ':'
          fldlist_lo_only(ipos+2:) = 'norm8wt'
       else
          fldlist_lo_only = 'norm8wt'
       end if
       nlo = nlo + 1
    end if

    ! Null-terminate for iMOAB.
    if (ncaas > 0) then
       ipos = len_trim(fldlist_caas)
       fldlist_caas(ipos+1:ipos+1) = C_NULL_CHAR
    else
       fldlist_caas = C_NULL_CHAR
    end if

    if (nlo > 0) then
       ipos = len_trim(fldlist_lo_only)
       fldlist_lo_only(ipos+1:ipos+1) = C_NULL_CHAR
    else
       fldlist_lo_only = C_NULL_CHAR
    end if
  end subroutine build_nlmap_sublists

  !=======================================================================

end module seq_map_mod
