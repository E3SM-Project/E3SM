module component_type_mod

  !----------------------------------------------------------------------------
  ! share code & libs
  !----------------------------------------------------------------------------
  use shr_kind_mod     , only: r8 => SHR_KIND_R8
  use shr_kind_mod     , only: cs => SHR_KIND_CS
  use shr_kind_mod     , only: cl => SHR_KIND_CL
  use shr_kind_mod     , only: IN => SHR_KIND_IN
  use seq_cdata_mod    , only: seq_cdata
  use seq_map_type_mod , only: seq_map
  use seq_comm_mct     , only: seq_comm_namelen
  use seq_comm_mct     , only: num_inst_atm, num_inst_lnd, num_inst_rof
  use seq_comm_mct     , only: num_inst_ocn, num_inst_ice, num_inst_glc
  use seq_comm_mct     , only: num_inst_wav, num_inst_esp, num_inst_iac
  use mct_mod

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  !
  ! on component pes
  public :: component_get_c2x_cc
  public :: component_get_x2c_cc
  public :: component_get_dom_cc
  public :: component_get_gsmap_cc
  public :: component_get_cdata_cc
  public :: component_get_iamroot_compid
  public :: check_fields
  !
  ! on cpl pes
  public :: component_get_x2c_cx
  public :: component_get_c2x_cx
  public :: component_get_dom_cx
  public :: component_get_gsmap_cx
  public :: component_get_drv2mdl
  public :: component_get_mdl2drv
  !
  ! on union coupler/component pes
  public :: component_get_mapper_Cc2x
  public :: component_get_mapper_Cx2c
  !
  ! on driver pes (all pes)
  public :: component_get_name
  public :: component_get_suffix
  public :: component_get_iamin_compid

! this is to replicate mct grid of a cx   
  public :: expose_mct_grid_moab
#ifdef MOABCOMP
  public :: compare_mct_av_moab_tag
#endif

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  type component_type
     !
     ! Coupler pes
     ! used by prep_xxx and all other coupler based routines
     !
     type(mct_ggrid) , pointer       :: dom_cx      => null() ! component domain (same for all instances)
     type(mct_gsMap) , pointer       :: gsMap_cx    => null() ! decomposition on coupler pes (same for all instances)
     type(mct_aVect) , pointer       :: x2c_cx      => null() !
     type(mct_aVect) , pointer       :: c2x_cx      => null()
     !
     ! Component pes
     !
     type(seq_cdata) , pointer       :: cdata_cc    => null()
     type(mct_ggrid) , pointer       :: dom_cc      => null()
     type(mct_gsMap) , pointer       :: gsMap_cc    => null() ! decomposition on component pes
     type(mct_aVect) , pointer       :: x2c_cc      => null()
     type(mct_aVect) , pointer       :: c2x_cc      => null()
     real(r8)        , pointer       :: drv2mdl(:)  => null() ! area correction factors
     real(r8)        , pointer       :: mdl2drv(:)  => null() ! area correction factors
#ifdef HAVE_MOAB
     integer                         :: mbApCCid ! moab app id on component side 
     integer                         :: mbGridType ! 0 for PC, 1 for cell (ocean, ice)  
     integer                         :: mblsize    ! size of local arrays
#endif 
     !
     ! Union of coupler/component pes - used by exchange routines
     !
     type(seq_map)   , pointer       :: mapper_Cc2x => null() ! coupler   -> component rearranging
     type(seq_map)   , pointer       :: mapper_Cx2c => null() ! component -> coupler   rearranging
     !
     ! Driver pes (all pes)
     !
     integer                         :: compid
     integer                         :: cplcompid
     integer                         :: cplallcompid
     integer                         :: mpicom_compid
     integer                         :: mpicom_cplcompid
     integer                         :: mpicom_cplallcompid
     logical                         :: iamin_compid
     logical                         :: iamin_cplcompid
     logical                         :: iamin_cplallcompid
     logical                         :: iamroot_compid
     logical                         :: present ! true => component is present and not stub
     integer                         :: nthreads_compid
     character(len=CL)               :: suffix
     character(len=1)                :: oneletterid
     character(len=3)                :: ntype
     character(len=seq_comm_namelen) :: name
  end type component_type

  public :: component_type

  !----------------------------------------------------------------------------
  ! Component type instances
  !----------------------------------------------------------------------------

  type(component_type), target :: atm(num_inst_atm)
  type(component_type), target :: lnd(num_inst_lnd)
  type(component_type), target :: rof(num_inst_rof)
  type(component_type), target :: ocn(num_inst_ocn)
  type(component_type), target :: ice(num_inst_ice)
  type(component_type), target :: glc(num_inst_glc)
  type(component_type), target :: wav(num_inst_wav)
  type(component_type), target :: esp(num_inst_esp)
  type(component_type), target :: iac(num_inst_iac)

  public :: atm, lnd, rof, ocn, ice, glc, wav, esp, iac

  !===============================================================================

contains

  !===============================================================================
  ! Accessor functions into component instance
  !===============================================================================

  function component_get_c2x_cc(comp)
    type(component_type), intent(in), target :: comp
    type(mct_avect), pointer    :: component_get_c2x_cc
    component_get_c2x_cc => comp%c2x_cc
  end function component_get_c2x_cc

  function component_get_c2x_cx(comp)
    type(component_type), intent(in), target :: comp
    type(mct_avect), pointer    :: component_get_c2x_cx
    component_get_c2x_cx => comp%c2x_cx
  end function component_get_c2x_cx

  function component_get_x2c_cc(comp)
    type(component_type), intent(in), target :: comp
    type(mct_avect), pointer    :: component_get_x2c_cc
    component_get_x2c_cc => comp%x2c_cc
  end function component_get_x2c_cc

  function component_get_x2c_cx(comp)
    type(component_type), intent(in), target :: comp
    type(mct_avect), pointer    :: component_get_x2c_cx
    component_get_x2c_cx => comp%x2c_cx
  end function component_get_x2c_cx

  function component_get_name(comp)
    type(component_type), intent(in), target :: comp
    character(len=seq_comm_namelen) :: component_get_name
    component_get_name = comp%name
  end function component_get_name

  function component_get_iamin_compid(comp)
    type(component_type), intent(in), target :: comp
    logical :: component_get_iamin_compid
    component_get_iamin_compid = comp%iamin_compid
  end function component_get_iamin_compid

  function component_get_iamroot_compid(comp)
    type(component_type), intent(in), target :: comp
    logical :: component_get_iamroot_compid
    component_get_iamroot_compid = comp%iamroot_compid
  end function component_get_iamroot_compid

  function component_get_suffix(comp)
    type(component_type), intent(in), target :: comp
    character(len=CL) :: component_get_suffix
    component_get_suffix = comp%suffix
  end function component_get_suffix

  function component_get_dom_cx(comp)
    type(component_type), intent(in), target :: comp
    type(mct_ggrid), pointer :: component_get_dom_cx
    component_get_dom_cx => comp%dom_cx
  end function component_get_dom_cx

  function component_get_dom_cc(comp)
    type(component_type), intent(in), target :: comp
    type(mct_ggrid), pointer :: component_get_dom_cc
    component_get_dom_cc => comp%dom_cc
  end function component_get_dom_cc

  function component_get_gsmap_cx(comp)
    type(component_type), intent(in), target :: comp
    type(mct_gsmap), pointer :: component_get_gsmap_cx
    component_get_gsmap_cx => comp%gsmap_cx
  end function component_get_gsmap_cx

  function component_get_gsmap_cc(comp)
    type(component_type), intent(in), target :: comp
    type(mct_gsmap), pointer :: component_get_gsmap_cc
    component_get_gsmap_cc => comp%gsmap_cc
  end function component_get_gsmap_cc

  function component_get_cdata_cc(comp)
    type(component_type), intent(in), target :: comp
    type(seq_cdata), pointer :: component_get_cdata_cc
    component_get_cdata_cc => comp%cdata_cc
  end function component_get_cdata_cc

  function component_get_drv2mdl(comp)
    type(component_type), intent(in), target :: comp
    real(r8), pointer :: component_get_drv2mdl(:)
    component_get_drv2mdl => comp%drv2mdl
  end function component_get_drv2mdl

  function component_get_mdl2drv(comp)
    type(component_type), intent(in), target :: comp
    real(r8), pointer :: component_get_mdl2drv(:)
    component_get_mdl2drv => comp%mdl2drv
  end function component_get_mdl2drv

  function component_get_mapper_Cc2x(comp)
    type(component_type), intent(in), target :: comp
    type(seq_map), pointer :: component_get_mapper_Cc2x
    component_get_mapper_Cc2x => comp%mapper_Cc2x
  end function component_get_mapper_Cc2x

  function component_get_mapper_Cx2c(comp)
    type(component_type), intent(in), target :: comp
    type(seq_map), pointer :: component_get_mapper_Cx2c
    component_get_mapper_Cx2c => comp%mapper_Cx2c
  end function component_get_mapper_Cx2c

  subroutine check_fields(comp, comp_index)
    use shr_infnan_mod, only: shr_infnan_isnan
    use mct_mod, only: mct_avect_getrlist2c, mct_gsMap_orderedPoints
    type(component_type), intent(in) :: comp
    integer(in), intent(in) :: comp_index

    integer(IN)   :: lsize             ! size of attr vect
    integer(IN)   :: nflds             ! number of attr vects
    integer(in)   :: fld, n            ! iterators
    integer(IN)  :: rank
    integer(IN) :: ierr
    integer(IN), pointer :: gpts(:)
    character(len=CL) :: msg

    if(associated(comp%c2x_cc) .and. associated(comp%c2x_cc%rattr)) then
       lsize = mct_avect_lsize(comp%c2x_cc)
       nflds = size(comp%c2x_cc%rattr,1)
       ! c2x_cc is allocated even if not used such as in stub models
       ! do not test this case.
       if(lsize <= 1 .and. nflds <= 1) return
#ifndef CPRFJ
       if(any(shr_infnan_isnan(comp%c2x_cc%rattr))) then
          do fld=1,nflds
             do n=1,lsize
                if(shr_infnan_isnan(comp%c2x_cc%rattr(fld,n))) then
                   call mpi_comm_rank(comp%mpicom_compid, rank, ierr)
                   call mct_gsMap_orderedPoints(comp%gsmap_cc, rank, gpts)
                   write(msg,'(a,a,a,i4,a,a,a,i8)')'component_mod:check_fields NaN found in ',trim(comp%name),' instance: ',&
                        comp_index,' field ',trim(mct_avect_getRList2c(fld, comp%c2x_cc)), ' 1d global index: ',gpts(n)
                   call shr_sys_abort(msg)
                endif
             enddo
          enddo
       endif
#endif
    endif
  end subroutine check_fields

  subroutine expose_mct_grid_moab (comp, imoabAPI)
    use shr_mpi_mod,       only: shr_mpi_commrank, shr_mpi_commsize
    use shr_sys_mod
    use shr_const_mod, only: SHR_CONST_PI
    use seq_comm_mct, only: cplid
    use seq_comm_mct, only: seq_comm_iamin
    use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs
    use iMOAB, only :  iMOAB_RegisterApplication, iMOAB_CreateVertices, iMOAB_WriteMesh, &
         iMOAB_DefineTagStorage, iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
         iMOAB_ResolveSharedEntities

    type(component_type), intent(in) :: comp
    integer , intent(out) :: imoabAPI

    integer                :: lsz
    type(mct_gGrid), pointer :: dom
    integer  :: mpicom_CPLID          ! MPI cpl communicator
    integer  :: iamcomp , iamcpl
    integer  :: ext_id

    ! local variables to fill in data
    integer, dimension(:), allocatable :: vgids
    !  retrieve everything we need from mct
    ! number of vertices is the size of mct grid
    real(r8), dimension(:), allocatable :: moab_vert_coords  ! temporary
    real(r8)   :: latv, lonv
    integer   dims, i, ilat, ilon, igdx, ierr, tagindex, ixarea, ixfrac
    integer tagtype, numco, ent_type
    character*100 outfile, wopts, localmeshfile, tagname
    character*32 appname

!----- formats -----
    character(*),parameter :: subName = '(expose_mct_grid_moab) '

    dims  = 3 ! store as 3d mesh

    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)
    
    if (seq_comm_iamin(CPLID)) then
      call shr_mpi_commrank(mpicom_CPLID, iamcpl  , 'expose_mct_grid_moab')
      dom => component_get_dom_cx(comp)
      lsz = mct_gGrid_lsize(dom)
      !print *, 'lsize: cx', lsz, ' iamcpl ' , iamcpl
      appname=comp%ntype//"CPMOAB"//CHAR(0)
      ! component instance
      ext_id = comp%compid + 200 ! avoid reuse
      ierr = iMOAB_RegisterApplication(appname, mpicom_CPLID, ext_id, imoabAPI)
      if (ierr > 0 )  &
         call shr_sys_abort(subname//'Error: cannot register moab app')
      allocate(moab_vert_coords(lsz*dims))
      allocate(vgids(lsz))
      ilat = MCT_GGrid_indexRA(dom,'lat')
      ilon = MCT_GGrid_indexRA(dom,'lon')
      igdx = MCT_GGrid_indexIA(dom,'GlobGridNum')
      do i = 1, lsz
        latv = dom%data%rAttr(ilat, i) *SHR_CONST_PI/180.
        lonv = dom%data%rAttr(ilon, i) *SHR_CONST_PI/180.
        moab_vert_coords(3*i-2)=COS(latv)*COS(lonv)
        moab_vert_coords(3*i-1)=COS(latv)*SIN(lonv)
        moab_vert_coords(3*i  )=SIN(latv)
        vgids(i) = dom%data%iAttr(igdx, i)
      enddo

      ierr = iMOAB_CreateVertices(imoabAPI, lsz*3, dims, moab_vert_coords)
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to create MOAB vertices in land model')

      tagtype = 0  ! dense, integer
      numco = 1
      tagname='GLOBAL_ID'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to retrieve GLOBAL_ID tag ')

      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage ( imoabAPI, tagname, lsz , ent_type, vgids)
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to set GLOBAL_ID tag ')

      ! MOAB TODO is this needed ? no vertices should be shared here, maybe just to set the part tag ?
      ierr = iMOAB_ResolveSharedEntities( imoabAPI, lsz, vgids );
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to resolve shared entities')

      !there are no shared entities, but we will set a special partition tag, in order to see the
      ! partitions ; it will be visible with a Pseudocolor plot in VisIt
      tagname='partition'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to create new partition tag ')

      vgids = iamcpl
      ierr = iMOAB_SetIntTagStorage ( imoabAPI, tagname, lsz , ent_type, vgids)
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to set partition tag ')

      ! use moab_vert_coords as a data holder for a frac tag and area tag that we will create
      !   on the vertices; do not allocate other data array
      !  do not be confused by this !
      ixfrac = MCT_GGrid_indexRA(dom,'frac')
      ixarea = MCT_GGrid_indexRA(dom,'area')
      tagname='frac'//CHAR(0)
      tagtype = 1 ! dense, double
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to create frac tag ')

      do i = 1, lsz
        moab_vert_coords(i) = dom%data%rAttr(ixfrac, i)
      enddo
      ierr = iMOAB_SetDoubleTagStorage ( imoabAPI, tagname, lsz , ent_type, moab_vert_coords)
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to set frac tag ')

      tagname='area'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to create area tag ')
      do i = 1, lsz
        moab_vert_coords(i) = dom%data%rAttr(ixarea, i) ! use the same doubles for second tag :)
      enddo

      ierr = iMOAB_SetDoubleTagStorage ( imoabAPI, tagname, lsz , ent_type, moab_vert_coords )
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to set area tag ')

      deallocate(moab_vert_coords)
      deallocate(vgids)
#ifdef MOABDEBUG
      !     write out the mesh file to disk, in parallel
      outfile = 'WHOLE_cx_'//comp%ntype//'.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(imoabAPI, outfile, wopts)
      if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to write the land mesh file')
#endif
    endif

  end subroutine expose_mct_grid_moab

#ifdef MOABCOMP
  ! assumes everything is on coupler pes here, to make sense
  subroutine compare_mct_av_moab_tag(comp, attrVect, mct_field, appId, tagname, ent_type, difference, first_time)
    
    use shr_mpi_mod,       only: shr_mpi_sum
    use shr_kind_mod,     only:  CXX => shr_kind_CXX
    use seq_comm_mct , only : CPLID, seq_comm_iamroot
    use seq_comm_mct, only:   seq_comm_setptrs
    use iMOAB, only : iMOAB_DefineTagStorage,  iMOAB_GetDoubleTagStorage, &
       iMOAB_SetDoubleTagStorageWithGid, iMOAB_GetMeshInfo
    
    use iso_c_binding 

    type(component_type), intent(in) :: comp
    integer , intent(in) :: appId, ent_type
    type(mct_aVect) , intent(in), pointer       :: attrVect
    character(*) , intent(in)       :: mct_field
    character(*) , intent(in)       :: tagname

    real(r8)      , intent(out)     :: difference
    logical , intent(in)            :: first_time

    real(r8)  :: differenceg ! global, reduced diff
    type(mct_ggrid), pointer    :: dom
    integer   :: kgg, mbSize, nloc, index_avfield

     ! moab
     integer                  :: tagtype, numco,  tagindex, ierr
     character(CXX)           :: tagname_mct
     integer ,    allocatable :: GlobalIds(:) ! used for setting values associated with ids
     
     real(r8) , allocatable :: values(:), mct_values(:)
     integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
     integer :: mpicom
     logical   :: iamroot


     character(*),parameter :: subName = '(compare_mct_av_moab_tag) '

     call seq_comm_setptrs(CPLID, mpicom=mpicom)

     nloc = mct_avect_lsize(attrVect)
     allocate(GlobalIds(nloc))
     allocate(values(nloc))
     dom => component_get_dom_cx(comp)
     kgg = mct_aVect_indexIA(dom%data ,"GlobGridNum" ,perrWith=subName)
     GlobalIds = dom%data%iAttr(kgg,:)

     index_avfield     = mct_aVect_indexRA(attrVect,trim(mct_field))
     values(:) = attrVect%rAttr(index_avfield,:) 

     tagname_mct = 'mct_'//trim(tagname)//C_NULL_CHAR

     
     tagtype = 1 ! dense, double
     numco = 1
     if (first_time) then
        ierr = iMOAB_DefineTagStorage(appId, tagname_mct, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
            call shr_sys_abort(subname//'Error: fail to define new tag for mct')
     endif 
     ierr = iMOAB_SetDoubleTagStorageWithGid ( appId, tagname_mct, nloc , ent_type, values, GlobalIds )
     if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to set new tags')

     deallocate(values)
     ! now start comparing tags after set
     ierr  = iMOAB_GetMeshInfo ( appId, nvert, nvise, nbl, nsurf, nvisBC );
     if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to get mesh info')
     if (ent_type .eq. 0) then
        mbSize = nvert(1)
     else if (ent_type .eq. 1) then
        mbSize = nvise(1)
     endif
     allocate(values(mbSize))
     allocate(mct_values(mbSize))

     ierr = iMOAB_GetDoubleTagStorage ( appId, tagname_mct, mbSize , ent_type, mct_values)
     if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to get mct tag values')
     ierr = iMOAB_GetDoubleTagStorage ( appId, tagname, mbSize , ent_type, values)
     if (ierr > 0 )  &
        call shr_sys_abort(subname//'Error: fail to get moab tag values')
      
     values  = mct_values - values

     difference = dot_product(values, values)
     differenceg = 0. ! intel complained; why ?
     call shr_mpi_sum(difference,differenceg,mpicom,subname)
     difference = sqrt(differenceg)
     iamroot = seq_comm_iamroot(CPLID)
     if ( iamroot ) then
        print * , subname, trim(comp%ntype), ' comp, difference on tag ', trim(tagname), ' = ', difference
        !call shr_sys_abort(subname//'differences between mct and moab values')
     endif
     deallocate(GlobalIds)
     deallocate(values)
     deallocate(mct_values)


  end subroutine compare_mct_av_moab_tag

#endif
end module component_type_mod
