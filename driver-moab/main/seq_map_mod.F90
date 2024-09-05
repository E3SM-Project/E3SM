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
  use mct_mod
  use seq_comm_mct
  use component_type_mod
  use seq_map_type_mod

  implicit none
  save
  private  ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: seq_map_init_rcfile     ! cpl pes
  public :: moab_map_init_rcfile    ! cpl pes
  public :: seq_map_init_rearrolap  ! cpl pes
  public :: seq_map_initvect        ! cpl pes
  public :: seq_map_map             ! cpl pes
  public :: seq_map_mapvect         ! cpl pes
  public :: seq_map_readdata        ! cpl pes

  interface seq_map_avNorm
     module procedure seq_map_avNormArr
     module procedure seq_map_avNormAvF
  end interface seq_map_avNorm

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
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="copy")

       if (mapid > 0 .and. .not. skip_match) then
          call seq_map_mappoint(mapid,mapper)
       else
          if(skip_match .and. seq_comm_iamroot(CPLID)) then
             write(logunit,'(A)') subname, 'skip_match true, force new map'
          endif
          call seq_map_mapinit(mapper,mpicom)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          mapper%gsmap_s => component_get_gsmap_cx(comp_s)
          mapper%gsmap_d => component_get_gsmap_cx(comp_d)
       endif

    elseif (samegrid) then
       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")

       if (mapid > 0 .and. .not. skip_match) then
          call seq_map_mappoint(mapid,mapper)
       else
          if(skip_match .and. seq_comm_iamroot(CPLID)) then
             write(logunit,'(A)') subname, 'skip_match true, force new map'
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

       call seq_map_mapmatch(mapid,gsMap_s=gsMap_s,gsMap_d=gsMap_d,mapfile=mapfile,strategy=maptype)

       if (mapid > 0 .and. .not. skip_match) then
          call seq_map_mappoint(mapid,mapper)
       else
          if(skip_match .and. seq_comm_iamroot(CPLID)) then
             write(logunit,'(A)') subname, 'skip_match true, force new map'
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


  subroutine moab_map_init_rcfile( mbappid, mbtsid, type_grid, comp_s, comp_d, &
    maprcfile, maprcname, maprctype, samegrid, string, esmf_map)

   use iMOAB, only: iMOAB_LoadMappingWeightsFromFile
   implicit none
   !-----------------------------------------------------
   !
   ! Arguments
   !
   type(integer)        ,intent(in)            :: mbappid  ! moab app id, identifing the map from source to target
   type(integer)        ,intent(in)            :: mbtsid   ! moab app id, identifying the target (now), for row based distribution
   type(integer)        ,intent(in)            :: type_grid ! 1 for SE, 2 for PC, 3 for FV; should be a member data
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
   !type(mct_gsmap), pointer    :: gsmap_s ! temporary pointers
   !type(mct_gsmap), pointer    :: gsmap_d ! temporary pointers
   integer(IN)                 :: mpicom
   character(CX)               :: mapfile
   character(CX)               :: mapfile_term
   character(CL)               :: maptype
   integer(IN)                 :: mapid
   character(CX)               :: sol_identifier !   /* "scalar", "flux", "custom" */
   integer                     :: ierr
   integer                     :: col_or_row ! 0 for row based, 1 for col based (we use row distribution now)


   character(len=*),parameter  :: subname = "(moab_map_init_rcfile) "
   !-----------------------------------------------------

   if (seq_comm_iamroot(CPLID) .and. present(string)) then
      write(logunit,'(A)') subname//' called for '//trim(string)
   endif

   call seq_comm_setptrs(CPLID, mpicom=mpicom)

   ! --- Initialize Smatp
   call shr_mct_queryConfigFile(mpicom,maprcfile,maprcname,mapfile,maprctype,maptype)
   !call shr_mct_sMatPInitnc(mapper%sMatp, mapper%gsMap_s, mapper%gsMap_d, trim(mapfile),trim(maptype),mpicom)
   sol_identifier = 'map-from-file'//CHAR(0)
   mapfile_term = trim(mapfile)//CHAR(0)
   if (seq_comm_iamroot(CPLID)) then
       write(logunit,*) subname,' reading map file with iMOAB: ', mapfile_term
   endif

   col_or_row = 0 ! row based distribution

   ierr = iMOAB_LoadMappingWeightsFromFile( mbappid, mbtsid, col_or_row, type_grid, sol_identifier, mapfile_term)
   if (ierr .ne. 0) then
      write(logunit,*) subname,' error in loading map file'
      call shr_sys_abort(subname//' ERROR in loading map file')
    endif
   if (seq_comm_iamroot(CPLID)) then
      write(logunit,'(2A,I6,4A)') subname,' iMOAB map app ID, maptype, mapfile = ', &
         mbappid,' ',trim(maptype),' ',trim(mapfile)
      call shr_sys_flush(logunit)
   endif

end subroutine moab_map_init_rcfile

  !=======================================================================

  subroutine seq_map_init_rearrolap(mapper, comp_s, comp_d, string, no_match)

    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(seq_map)        ,intent(inout),pointer :: mapper
    type(component_type) ,intent(inout)         :: comp_s
    type(component_type) ,intent(inout)         :: comp_d
    character(len=*)     ,intent(in),optional   :: string
    logical              ,intent(in),optional   :: no_match
    !
    ! Local Variables
    !
    integer(IN)                :: mapid
    type(mct_gsmap), pointer   :: gsmap_s
    type(mct_gsmap), pointer   :: gsmap_d
    integer(IN)                :: mpicom
    logical                    :: skip_match
    character(len=*),parameter :: subname = "(seq_map_init_rearrolap) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    call seq_comm_setptrs(CPLID, mpicom=mpicom)

    gsmap_s => component_get_gsmap_cx(comp_s)
    gsmap_d => component_get_gsmap_cx(comp_d)

    skip_match = .false.
    if (present(no_match)) then
       if (no_match) skip_match = .true.
    endif
    if (mct_gsmap_Identical(gsmap_s,gsmap_d) .and. .not.skip_match ) then
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

    use iso_c_binding
    use iMOAB, only: iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, iMOAB_SetDoubleTagStorage, &
      iMOAB_GetIntTagStorage, iMOAB_ApplyScalarProjectionWeights, &
      iMOAB_SendElementTag, iMOAB_ReceiveElementTag, iMOAB_FreeSenderBuffers

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
#ifdef HAVE_MOAB
    logical  :: valid_moab_context
    integer  :: ierr, nfields, lsize_src, lsize_tgt, arrsize_tgt, j, arrsize_src
    character(len=CXX) :: fldlist_moab
    character(len=CXX) :: tagname
    integer    :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
    type(mct_list) :: temp_list
    integer, dimension(:), allocatable  :: globalIds
    real(r8), dimension(:), allocatable  :: wghts
    real(kind=r8) , allocatable  :: targtags(:,:), targtags_ini(:,:)
    real(kind=r8)  :: factor
    integer  :: filter_type ! used for caas projection
#endif
    !
    ! Local Variables
    !
    logical :: lnorm  ! true if normalization is to be done
    logical :: mbnorm ! moab copy of lnorm
    logical :: mbpresent ! moab logical for presence of norm weight string
    integer(IN),save :: ltag    ! message tag for rearrange
    character(len=*),parameter :: subname = "(seq_map_map) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

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

#ifdef HAVE_MOAB
       ! check whether the application ID is defined on the current process
       if ( mapper%src_mbid .lt. 0 .or. mapper%tgt_mbid .lt. 0 ) then
         valid_moab_context = .FALSE.
       else
         valid_moab_context = .TRUE.
       endif

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

         if (mbnorm) then
           fldlist_moab = trim(fldlist_moab)//":norm8wt"//C_NULL_CHAR
           nfields=nfields + 1
         else
           fldlist_moab = trim(fldlist_moab)//C_NULL_CHAR
         endif


#ifdef MOABDEBUG
         if (seq_comm_iamroot(CPLID)) then
            write(logunit,*) subname, 'iMOAB mapper ',trim(mapper%mbname), ' iMOAB_mapper  nfields', &
                  nfields,  ' fldlist_moab=', trim(fldlist_moab)
            call shr_sys_flush(logunit)
         endif
#endif
       endif ! valid_moab_context
#endif

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

    endif

    if (mapper%copy_only .or. mapper%rearrange_only) then

#ifdef HAVE_MOAB
       if ( valid_moab_context ) then
#ifdef MOABDEBUG
         if (seq_comm_iamroot(CPLID)) then
            write(logunit, *) subname,' iMOAB mapper rearrange or copy ', mapper%mbname, ' send/recv tags ', trim(fldlist_moab), &
              ' mbpresent=', mbpresent, ' mbnorm=', mbnorm
            call shr_sys_flush(logunit)
         endif
#endif
         ierr = iMOAB_SendElementTag( mapper%src_mbid, fldlist_moab, mapper%mpicom, mapper%intx_context );
         if (ierr .ne. 0) then
            write(logunit, *) subname,' iMOAB mapper ', mapper%mbname, ' error in sending tags ', trim(fldlist_moab), ierr
            call shr_sys_flush(logunit)
            call shr_sys_abort(subname//' ERROR in sending tags')
         endif
         ! receive in the target app
         ierr = iMOAB_ReceiveElementTag( mapper%tgt_mbid, fldlist_moab, mapper%mpicom, mapper%src_context );
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in receiving tags iMOAB mapper ', mapper%mbname,  trim(fldlist_moab)
            call shr_sys_flush(logunit)
            call shr_sys_abort(subname//' ERROR in receiving tags')
         endif
         ! now free buffers
         ierr = iMOAB_FreeSenderBuffers( mapper%src_mbid, mapper%intx_context )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in freeing buffers ', trim(fldlist_moab)
            call shr_sys_abort(subname//' ERROR in freeing buffers') ! serious enough
         endif
       endif ! if (valid_moab_context)

#endif

      else

#ifdef HAVE_MOAB
       if ( valid_moab_context ) then
       ! NORMALIZATION
         if (mbnorm .or. mbpresent) then
            !  initialize the weight tag and multiply it by the input tags.
            ! get target mesh info
            ierr  = iMOAB_GetMeshInfo ( mapper%src_mbid, nvert, nvise, nbl, nsurf, nvisBC );
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error getting mesh info for ', mapper%mbname
               call shr_sys_abort(subname//' ERROR getting mesh info') ! serious enough
            endif
            lsize_src = nvise(1) ! number of active cells

            ! init normalization weight
            allocate(wghts(lsize_src))
            wghts = 1.0_r8
            tagname = "norm8wt"//C_NULL_CHAR
            ! set the normalization factor to 1
            ierr = iMOAB_SetDoubleTagStorage (mapper%src_mbid, tagname, lsize_src , mapper%tag_entity_type, wghts)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error setting init value for mapping norm factor ',ierr,trim(tagname)
               call shr_sys_abort(subname//' ERROR setting norm init value') ! serious enough
            endif
#ifdef MOABDEBUG
            if (seq_comm_iamroot(CPLID)) then
               write(logunit, *) subname,' iMOAB mapper ', mapper%mbname, ' set norm8wt 1  on source with app id: ', mapper%src_mbid
               call shr_sys_flush(logunit)
            endif
#endif
            ! if a normalization factor was specified, get it and multiply src tags by it
            if(mbpresent) then
               tagname = avwtsfld_s//C_NULL_CHAR
               ierr = iMOAB_GetDoubleTagStorage (mapper%src_mbid, tagname, lsize_src , mapper%tag_entity_type, wghts)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error getting value for mapping norm factor ', trim(tagname)
                  call shr_sys_abort(subname//' ERROR getting norm factor') ! serious enough
               endif

               ! get the fieldlist including weight
               allocate(targtags(lsize_src,nfields))
               allocate(targtags_ini(lsize_src,nfields))
               arrsize_src=lsize_src*(nfields)

               ! get the current values of all source tags including the norm8wt currently set to 1
               ierr = iMOAB_GetDoubleTagStorage (mapper%src_mbid, fldlist_moab, arrsize_src , mapper%tag_entity_type, targtags)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error getting source tag values ', mapper%mbname, mapper%src_mbid, trim(fldlist_moab), arrsize_src, mapper%tag_entity_type
                  call shr_sys_abort(subname//' ERROR getting source tag values') ! serious enough
               endif

               targtags_ini = targtags
               ! multiply by the value of the avwtsfld_s field.
               ! norm8wt is 1 so it will record the value of the weight.
               do j = 1, lsize_src
                 targtags(j,:)= targtags(j,:)*wghts(j)
               enddo
#ifdef MOABDEBUG
         if (seq_comm_iamroot(CPLID)) then
            write(logunit, *) subname,' iMOAB projection mapper: ', mapper%mbname, ' normalize nfields=', &
               nfields, ' arrsize_src on root:', arrsize_src, ' shape(targtags_ini)=', shape(targtags_ini)
            call shr_sys_flush(logunit)
         endif
#endif
               ! put the new values on the mesh for later mapping
               ierr = iMOAB_SetDoubleTagStorage (mapper%src_mbid, fldlist_moab, arrsize_src , mapper%tag_entity_type, targtags)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error setting normed source tag values ', mapper%mbname
                  call shr_sys_abort(subname//' ERROR setting normed source tag values') ! serious enough
               endif
               deallocate(targtags)
            endif ! end multiplication by norm factor
            deallocate(wghts)
         endif  ! end NORMALIZATION

         !
         ierr = iMOAB_SendElementTag( mapper%src_mbid, fldlist_moab, mapper%mpicom, mapper%intx_context );
         if (ierr .ne. 0) then
            write(logunit, *) subname,' iMOAB mapper error in sending tags ', mapper%mbname,  trim(fldlist_moab)
            call shr_sys_flush(logunit)
            call shr_sys_abort(subname//' ERROR in sending tags')
         endif
       endif
       if ( valid_moab_context ) then
         ! receive in the intx app, because it is redistributed according to coverage (trick)
         ! for true intx cases, tgt_mbid is set to be the same as intx_mbid
         ! just read map is special
         if (mapper%read_map)  then ! receive indeed in target app
#ifdef MOABDEBUG
            if (seq_comm_iamroot(CPLID)) then
               write(logunit, *) subname,' iMOAB mapper receiving tags with read_map and tgt_mbid: ', &
                mapper%mbname, trim(fldlist_moab)
            endif
#endif
            ierr = iMOAB_ReceiveElementTag( mapper%tgt_mbid, fldlist_moab, mapper%mpicom, mapper%src_context )
         else ! receive in the intx app, trick
#ifdef MOABDEBUG
            if (seq_comm_iamroot(CPLID)) then
               write(logunit, *) subname,' iMOAB mapper receiving tags with intx and intx_mbid: ', &
                mapper%mbname, trim(fldlist_moab)
            endif
#endif
            ierr = iMOAB_ReceiveElementTag( mapper%intx_mbid, fldlist_moab, mapper%mpicom, mapper%src_context )
         endif
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in receiving tags ', mapper%mbname, 'recv:',  mapper%intx_mbid, trim(fldlist_moab)
            call shr_sys_flush(logunit)
            call shr_sys_abort(subname//' ERROR in receiving tags')
            !valid_moab_context = .false. ! do not attempt to project
         endif
         ! now free buffers
         ierr = iMOAB_FreeSenderBuffers( mapper%src_mbid, mapper%intx_context )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in freeing buffers ', trim(fldlist_moab)
            call shr_sys_abort(subname//' ERROR in freeing buffers') ! serious enough
         endif
       endif
       if ( valid_moab_context ) then

#ifdef MOABDEBUG
         if (seq_comm_iamroot(CPLID)) then
            write(logunit, *) subname,' iMOAB projection mapper: ',trim(mapper%mbname), ' between ', mapper%src_mbid, ' and ',  mapper%tgt_mbid, trim(fldlist_moab)
            call shr_sys_flush(logunit)
         endif
#endif
         filter_type = 0 ! no
         ierr = iMOAB_ApplyScalarProjectionWeights ( mapper%intx_mbid, filter_type, mapper%weight_identifier, fldlist_moab, fldlist_moab)
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in applying weights '
            call shr_sys_abort(subname//' ERROR in applying weights')
         endif

         ! complete the normalization process
         if (mbnorm) then
            ierr  = iMOAB_GetMeshInfo ( mapper%tgt_mbid, nvert, nvise, nbl, nsurf, nvisBC );
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error getting mesh info for target ', mapper%mbname
               call shr_sys_abort(subname//' ERROR getting mesh info') ! serious enough
            endif

            lsize_tgt = nvise(1) ! number of active cells
            tagname = "norm8wt"//C_NULL_CHAR
            allocate(wghts(lsize_tgt))

            ! get values of weights after mapping
            ierr = iMOAB_GetDoubleTagStorage (mapper%tgt_mbid, tagname, lsize_tgt , mapper%tag_entity_type, wghts)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error getting value for mapping norm factor post-map ', ierr, trim(tagname)
               call shr_sys_abort(subname//' ERROR getting norm factor') ! serious enough
            endif

            ! get values of target tags after mapping
            allocate(targtags(lsize_tgt,nfields))
            arrsize_tgt=lsize_tgt*(nfields)
            ierr = iMOAB_GetDoubleTagStorage (mapper%tgt_mbid, fldlist_moab, arrsize_tgt , mapper%tag_entity_type, targtags)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error getting destination tag values ', mapper%mbname
               call shr_sys_abort(subname//' ERROR getting source tag values') ! serious enough
            endif

            ! do the post mapping normalization
            ! TODO:  add some check for wghts < puny
            do j = 1, lsize_tgt
               factor = wghts(j)
               if (wghts(j) .ne. 0) factor = 1.0_r8/wghts(j) ! should we compare to a small value instead ?
               targtags(j,:)= targtags(j,:)*factor
            enddo

            ! put the values back on the mesh
            ierr = iMOAB_SetDoubleTagStorage (mapper%tgt_mbid, fldlist_moab, arrsize_tgt , mapper%tag_entity_type, targtags)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error getting destination tag values ', mapper%mbname
               call shr_sys_abort(subname//' ERROR getting source tag values') ! serious enough
            endif

            deallocate(wghts, targtags)
            if (mbpresent) then
#ifdef MOABDEBUG
               if (seq_comm_iamroot(CPLID)) then
                  write(logunit, *) subname,' iMOAB projection mapper: ', mapper%mbname, ' shape(targtags_ini)=', shape(targtags_ini)
                  call shr_sys_flush(logunit)
               endif
#endif
               ! put the values back on the source mesh
               ierr = iMOAB_SetDoubleTagStorage (mapper%src_mbid, fldlist_moab, arrsize_src , mapper%tag_entity_type, targtags_ini)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error setting source tag values ', mapper%mbname
                  call shr_sys_abort(subname//' ERROR setting source tag values') ! serious enough
               endif
               deallocate(targtags_ini)
            endif

         endif ! end normalization

       endif
#endif

      endif ! end of mapping type if else

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
    character(len=CL)          :: lstring
    character(len=*),parameter :: subname = "(seq_map_initvect) "
    !-----------------------------------------------------

    lstring = ' '
    if (present(string)) then
       if (seq_comm_iamroot(CPLID)) write(logunit,'(A)') subname//' called for '//trim(string)
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
       write(logunit,'(A)') subname//' called for '//trim(string)
       call shr_sys_flush(logunit)
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
    integer(IN) :: lsize_i, lsize_f, kf, j
    real(r8),allocatable :: frac_i(:)
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
    type(mct_aVect)        :: avp_i , avp_o
    integer(IN)            :: j,kf
    integer(IN)            :: lsize_i,lsize_o
    real(r8)               :: normval
    character(CX)          :: lrList,appnd
    logical                :: lnorm
    character(*),parameter :: subName = '(seq_map_avNormArr) '
    character(len=*),parameter :: ffld = 'norm8wt'  ! want something unique
    !-----------------------------------------------------

    lsize_i = mct_aVect_lsize(av_i)
    lsize_o = mct_aVect_lsize(av_o)

    lnorm = .true.
    if (present(norm)) then
       lnorm = norm
    endif

    if (present(norm_i)) then
       if (.not.lnorm) call shr_sys_abort(subname//' ERROR norm_i and norm = false')
       if (size(norm_i) /= lsize_i) call shr_sys_abort(subname//' ERROR size(norm_i) ne lsize_i')
    endif

    !--- create temporary avs for mapping ---

    if (lnorm .or. present(norm_i)) then
       appnd = ':'//ffld
    else
       appnd = ''
    endif
    if (present(rList)) then
       call mct_aVect_init(avp_i, rList=trim( rList)//trim(appnd), lsize=lsize_i)
       call mct_aVect_init(avp_o, rList=trim( rList)//trim(appnd), lsize=lsize_o)
    else
       lrList = ''
       if(mct_aVect_nRAttr(av_i) /= 0) lrList = mct_aVect_exportRList2c(av_i)
       call mct_aVect_init(avp_i, rList=trim(lrList)//trim(appnd), lsize=lsize_i)
       lrList = ''
       if(mct_aVect_nRAttr(av_o) /= 0) lrList = mct_aVect_exportRList2c(av_o)
       call mct_aVect_init(avp_o, rList=trim(lrList)//trim(appnd), lsize=lsize_o)
    endif

    !--- copy av_i to avp_i and set ffld value to 1.0
    !--- then multiply all fields by norm_i if norm_i exists
    !--- this will do the right thing for the norm_i normalization

    call mct_aVect_copy(aVin=av_i, aVout=avp_i, VECTOR=mct_usevector)
    if (lnorm .or. present(norm_i)) then
       kf = mct_aVect_indexRA(avp_i,ffld)
       do j = 1,lsize_i
          avp_i%rAttr(kf,j) = 1.0_r8
       enddo

       if (present(norm_i)) then
          !$omp simd
          do j = 1,lsize_i
             avp_i%rAttr(:,j) = avp_i%rAttr(:,j)*norm_i(j)
          enddo
       endif
    endif

    !--- map ---

    if (mapper%esmf_map) then
       call shr_sys_abort(subname//' ERROR: esmf SMM not supported')
    else
       ! MCT based SMM
       call mct_sMat_avMult(avp_i, mapper%sMatp, avp_o, VECTOR=mct_usevector)
    endif

    !--- renormalize avp_o by mapped norm_i  ---

    if (lnorm) then
       kf = mct_aVect_indexRA(avp_o,ffld)
       !$omp simd
       do j = 1,lsize_o
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
