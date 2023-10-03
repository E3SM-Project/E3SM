module seq_map_type_mod

  use shr_kind_mod , only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod , only: CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use shr_mct_mod  , only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use shr_sys_mod
  use shr_const_mod
  use seq_comm_mct,  only: logunit, CPLID, seq_comm_iamroot
  use mct_mod

  type seq_map
     logical                 :: copy_only
     logical                 :: rearrange_only
     logical                 :: esmf_map
     type(mct_rearr)         :: rearr
     type(mct_sMatp)         :: sMatp
     !
     !---- for comparing
     integer(IN)             :: counter   ! indicates which seq_maps this mapper points to
     character(CL)           :: strategy  ! indicates the strategy for this mapper, (copy, rearrange, X, Y)
     character(CX)           :: mapfile   ! indicates the mapping file used
     type(mct_gsMap),pointer :: gsmap_s
     type(mct_gsMap),pointer :: gsmap_d
     !
     !---- for cart3d
     character(CL)           :: cart3d_init
     real(R8), pointer       :: slon_s(:)
     real(R8), pointer       :: clon_s(:)
     real(R8), pointer       :: slat_s(:)
     real(R8), pointer       :: clat_s(:)
     real(R8), pointer       :: slon_d(:)
     real(R8), pointer       :: clon_d(:)
     real(R8), pointer       :: slat_d(:)
     real(R8), pointer       :: clat_d(:)
     integer(IN)             :: mpicom    ! mpicom

#ifdef HAVE_MOAB
     ! MOAB additional members, that store source, target and intx MOAB appids
     !  these are integers greater than or equal to 0
     !
     !   in general, rearrange can be solved with a parcommgraph, true fv-fv intx with an actual intx
     !   and a weight matrix;
     ! intx appid should be used for one map usually
     ! source and target app ids also make sense only on the coupler pes
     integer                 :: src_mbid, tgt_mbid, intx_mbid, src_context, intx_context
     character*32            :: weight_identifier ! 'state' OR 'flux'
     character*16            :: mbname
     integer                 :: tag_entity_type
     integer                 :: nentities ! this should be used only if copy_only is true
     logical                 :: read_map
     !
#endif

  end type seq_map
  public seq_map

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! seq_map_maxcnt is the total number of mappings supported
  ! seq_map_cnt is the total number of mappings initialized any any time
  ! seq_maps are the mappers that have been initialized

  integer(IN),parameter :: seq_map_maxcnt = 5000
  integer(IN)           :: seq_map_cnt = 0
  type(seq_map),private,target  :: seq_maps(seq_map_maxcnt)

  !  tcraig, work-in-progress
  !  type seq_map_node
  !     type(seq_map_node), pointer :: next,prev
  !     type(seq_map), pointer :: seq_map
  !  end type seq_map_node
  !  type(seq_map_node), pointer :: seq_map_list, seq_map_curr

  !===============================================================================
contains
  !===============================================================================

  subroutine seq_map_mapmatch(mapid,gsMap_s,gsMap_d,mapfile,strategy)

    ! This method searches through the current seq_maps to find a
    ! mapping file that matches the values passed in

    implicit none
    integer         ,intent(out) :: mapid
    type(mct_gsMap) ,intent(in),optional :: gsMap_s
    type(mct_gsMap) ,intent(in),optional :: gsMap_d
    character(len=*),intent(in),optional :: mapfile
    character(len=*),intent(in),optional :: strategy

    integer(IN) :: m
    logical     :: match
    character(*),parameter :: subName = '(seq_map_mapmatch) '

    mapid = -1
    ! tcraig - this return turns off the mapping reuse
    ! RETURN

    do m = 1,seq_map_cnt
       match = .true.

       if (match .and. present(mapfile)) then
          if (trim(mapfile) /= trim(seq_maps(m)%mapfile)) match = .false.
       endif
       if (match .and. present(strategy)) then
          if (trim(strategy) /= trim(seq_maps(m)%strategy)) match = .false.
       endif
       if (match .and. present(gsMap_s)) then
          if (.not.mct_gsmap_Identical(gsmap_s,seq_maps(m)%gsmap_s)) match = .false.
       endif
       if (match .and. present(gsMap_d)) then
          if (.not.mct_gsmap_Identical(gsmap_d,seq_maps(m)%gsmap_d)) match = .false.
       endif

       if (match) then
          mapid = m
          if (seq_comm_iamroot(CPLID)) then
             write(logunit,'(A,i6)') subname//' found match ',mapid
#ifdef MOABCOMP
             write(logunit,'(A,i6)') subname//' strategy '//trim(seq_maps(mapid)%strategy)//&
             ' mapfile '//trim(seq_maps(mapid)%mapfile)
#endif
             call shr_sys_flush(logunit)
          endif
          return
       endif
    enddo

  end subroutine seq_map_mapmatch

  !===============================================================================

  subroutine seq_map_mapinit(mapper,mpicom)

    ! This method initializes a new seq_maps map datatype and
    ! has the mapper passed in point to it

    implicit none
    type(seq_map)   ,intent(inout),pointer           :: mapper
    integer(IN)     ,intent(in)                      :: mpicom

    character(*),parameter :: subName = '(seq_map_mapinit) '

    ! set the seq_map data
    seq_map_cnt = seq_map_cnt + 1
    if (seq_map_cnt > seq_map_maxcnt) then
       write(logunit,*) trim(subname),'seq_map_cnt too large',seq_map_cnt
       call shr_sys_abort(subName // "seq_map_cnt bigger than seq_map_maxcnt")
    endif
    mapper => seq_maps(seq_map_cnt)
    mapper%counter = seq_map_cnt

    mapper%copy_only      = .false.
    mapper%rearrange_only = .false.
    mapper%mpicom         = mpicom
    mapper%strategy       = "undefined"
    mapper%mapfile        = "undefined"

#ifdef HAVE_MOAB
    mapper%src_mbid  = -1
    mapper%tgt_mbid  = -1
    mapper%intx_mbid = -1
    mapper%nentities =  0
    mapper%tag_entity_type = 1 ! cells most of the time when we need it
    mapper%mbname    = "undefined"
    mapper%read_map = .false.
#ifdef MOABCOMP
    if (seq_comm_iamroot(CPLID)) then
      write(logunit,'(A,i6)') subname//' call init map for mapper with id ',mapper%counter
      call shr_sys_flush(logunit)
    endif
#endif
#endif

  end subroutine seq_map_mapinit

  !===============================================================================

  subroutine seq_map_mappoint(mapid,mapper)

    ! This method searches through the current seq_maps to find a
    ! mapping file that matches the values passed in

    implicit none
    integer         ,intent(in) :: mapid
    type(seq_map)   ,intent(inout),pointer :: mapper

    mapper => seq_maps(mapid)

  end subroutine seq_map_mappoint

  !===============================================================================

  subroutine seq_map_gsmapcheck(gsmap1,gsmap2)

    ! This method verifies that two gsmaps are of the same global size

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

  subroutine seq_map_set_type(mapper, mbid, ent_type)
    use iMOAB, only: iMOAB_GetMeshInfo
    type(seq_map)   ,intent(in),pointer :: mapper
    integer         ,intent(in) :: mbid
    integer         ,intent(in) :: ent_type
    
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3), ierr


    ierr  = iMOAB_GetMeshInfo ( mbid, nvert, nvise, nbl, nsurf, nvisBC ); 
    if (ent_type .eq. 0) then
       mapper%nentities = nvert(1)
    else if (ent_type .eq. 1) then
       mapper%nentities = nvise(1)
    endif

    mapper%tag_entity_type = ent_type

  end  subroutine seq_map_set_type

end module seq_map_type_mod
