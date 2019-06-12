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
  use seq_comm_mct     , only: num_inst_wav, num_inst_esp
  use mct_mod
  use seq_comm_mct     , only: CPLID
  use seq_comm_mct     , only: seq_comm_getinfo => seq_comm_setptrs
  use abortutils       , only : endrun
  implicit none
  save
  private
#include <mpif.h>

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
#ifdef MOABDEBUGMCT
  public :: expose_mct_grid_moab
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
     integer                         :: instn
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

  public :: atm, lnd, rof, ocn, ice, glc, wav, esp

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
    endif
  end subroutine check_fields

#ifdef MOABDEBUGMCT
  subroutine expose_mct_grid_moab (comp)
    use shr_mpi_mod,       only: shr_mpi_commrank, shr_mpi_commsize
    type(component_type), intent(in) :: comp
    integer                :: lsz
    type(mct_gGrid), pointer :: dom
    integer  :: mpicom_CPLID          ! MPI cpl communicator
    integer  :: imoabAPI
    integer  :: iamcomp , iamcpl
    integer  :: ext_id
    integer , external :: iMOAB_RegisterFortranApplication, iMOAB_CreateVertices, iMOAB_WriteMesh, &
         iMOAB_DefineTagStorage, iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
         iMOAB_ResolveSharedEntities
    ! local variables to fill in data
    integer, dimension(:), allocatable :: vgids
    !  retrieve everything we need from land domain mct_ldom
    ! number of vertices is the size of land domain
    real(r8), dimension(:), allocatable :: moab_vert_coords  ! temporary
    real(r8)   :: latv, lonv
    integer   dims, i, ilat, ilon, igdx, ierr, tagindex, ixarea, ixfrac
    integer tagtype, numco, ent_type
    character*100 outfile, wopts, localmeshfile, tagname
    character*32 appname
    real(R8),parameter :: SHR_CONST_PI      = 3.14159265358979323846_R8  ! pi

    dims  =3 ! store as 3d mesh


    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)
    if (comp%iamin_compid) then
      call shr_mpi_commrank(comp%mpicom_compid, iamcomp  , 'expose_mct_grid_moab')
      dom => component_get_dom_cc(comp)
      lsz = mct_gGrid_lsize(dom)
      !print *, 'lsize: cc', lsz, ' iamcomp ' ,iamcomp
      appname=comp%ntype//"MOAB"//CHAR(0)
      ! component instance
      ext_id = comp%compid + 100 ! avoid reuse
      ierr = iMOAB_RegisterFortranApplication(appname, comp%mpicom_compid, ext_id, imoabAPI)
      if (ierr > 0 )  &
         call endrun('Error: cannot register moab app')
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
        call endrun('Error: fail to create MOAB vertices in land model')

      tagtype = 0  ! dense, integer
      numco = 1
      tagname='GLOBAL_ID'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to retrieve GLOBAL_ID tag ')

      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage ( imoabAPI, tagname, lsz , ent_type, vgids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set GLOBAL_ID tag ')

      ierr = iMOAB_ResolveSharedEntities( imoabAPI, lsz, vgids );
      if (ierr > 0 )  &
        call endrun('Error: fail to resolve shared entities')

      !there are no shared entities, but we will set a special partition tag, in order to see the
      ! partitions ; it will be visible with a Pseudocolor plot in VisIt
      tagname='partition'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create new partition tag ')

      vgids = iamcomp
      ierr = iMOAB_SetIntTagStorage ( imoabAPI, tagname, lsz , ent_type, vgids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set partition tag ')

      ! use moab_vert_coords as a data holder for a frac tag and area tag that we will create
      !   on the vertices; do not allocate other data array
      !  do not be confused by this !
      ixfrac = MCT_GGrid_indexRA(dom,'frac')
      ixarea = MCT_GGrid_indexRA(dom,'area')
      tagname='frac'//CHAR(0)
      tagtype = 1 ! dense, double
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create frac tag ')

      do i = 1, lsz
        moab_vert_coords(i) = dom%data%rAttr(ixfrac, i)
      enddo
      ierr = iMOAB_SetDoubleTagStorage ( imoabAPI, tagname, lsz , ent_type, moab_vert_coords)
      if (ierr > 0 )  &
        call endrun('Error: fail to set frac tag ')

      tagname='area'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create area tag ')
      do i = 1, lsz
        moab_vert_coords(i) = dom%data%rAttr(ixarea, i) ! use the same doubles for second tag :)
      enddo

      ierr = iMOAB_SetDoubleTagStorage ( imoabAPI, tagname, lsz , ent_type, moab_vert_coords )
      if (ierr > 0 )  &
        call endrun('Error: fail to set area tag ')

      deallocate(moab_vert_coords)
      deallocate(vgids)
      !     write out the mesh file to disk, in parallel
      outfile = 'WHOLE_'//comp%ntype//'.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(imoabAPI, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the land mesh file')
    endif
    if (mpicom_CPLID /= MPI_COMM_NULL) then
      call shr_mpi_commrank(mpicom_CPLID, iamcpl  , 'expose_mct_grid_moab')
      dom => component_get_dom_cx(comp)
      lsz = mct_gGrid_lsize(dom)
      !print *, 'lsize: cx', lsz, ' iamcpl ' , iamcpl
      appname=comp%ntype//"CPMOAB"//CHAR(0)
      ! component instance
      ext_id = comp%compid + 200 ! avoid reuse
      ierr = iMOAB_RegisterFortranApplication(appname, mpicom_CPLID, ext_id, imoabAPI)
      if (ierr > 0 )  &
         call endrun('Error: cannot register moab app')
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
        call endrun('Error: fail to create MOAB vertices in land model')

      tagtype = 0  ! dense, integer
      numco = 1
      tagname='GLOBAL_ID'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to retrieve GLOBAL_ID tag ')

      ent_type = 0 ! vertex type
      ierr = iMOAB_SetIntTagStorage ( imoabAPI, tagname, lsz , ent_type, vgids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set GLOBAL_ID tag ')

      ierr = iMOAB_ResolveSharedEntities( imoabAPI, lsz, vgids );
      if (ierr > 0 )  &
        call endrun('Error: fail to resolve shared entities')

      !there are no shared entities, but we will set a special partition tag, in order to see the
      ! partitions ; it will be visible with a Pseudocolor plot in VisIt
      tagname='partition'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create new partition tag ')

      vgids = iamcpl
      ierr = iMOAB_SetIntTagStorage ( imoabAPI, tagname, lsz , ent_type, vgids)
      if (ierr > 0 )  &
        call endrun('Error: fail to set partition tag ')

      ! use moab_vert_coords as a data holder for a frac tag and area tag that we will create
      !   on the vertices; do not allocate other data array
      !  do not be confused by this !
      ixfrac = MCT_GGrid_indexRA(dom,'frac')
      ixarea = MCT_GGrid_indexRA(dom,'area')
      tagname='frac'//CHAR(0)
      tagtype = 1 ! dense, double
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create frac tag ')

      do i = 1, lsz
        moab_vert_coords(i) = dom%data%rAttr(ixfrac, i)
      enddo
      ierr = iMOAB_SetDoubleTagStorage ( imoabAPI, tagname, lsz , ent_type, moab_vert_coords)
      if (ierr > 0 )  &
        call endrun('Error: fail to set frac tag ')

      tagname='area'//CHAR(0)
      ierr = iMOAB_DefineTagStorage(imoabAPI, tagname, tagtype, numco,  tagindex )
      if (ierr > 0 )  &
        call endrun('Error: fail to create area tag ')
      do i = 1, lsz
        moab_vert_coords(i) = dom%data%rAttr(ixarea, i) ! use the same doubles for second tag :)
      enddo

      ierr = iMOAB_SetDoubleTagStorage ( imoabAPI, tagname, lsz , ent_type, moab_vert_coords )
      if (ierr > 0 )  &
        call endrun('Error: fail to set area tag ')

      deallocate(moab_vert_coords)
      deallocate(vgids)
      !     write out the mesh file to disk, in parallel
      outfile = 'WHOLE_cx_'//comp%ntype//'.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(imoabAPI, outfile, wopts)
      if (ierr > 0 )  &
        call endrun('Error: fail to write the land mesh file')
    endif

  end subroutine expose_mct_grid_moab
#endif
end module component_type_mod
