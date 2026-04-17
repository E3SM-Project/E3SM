module cplcomp_exchange_mod

  use shr_kind_mod, only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod, only: CL => SHR_KIND_CL, CX => SHR_KIND_CX, CXX => SHR_KIND_CXX
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod,  only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use mct_mod
  use seq_map_type_mod
  use component_type_mod
  use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other
  use seq_flds_mod, only: seq_flds_dom_fields
  use seq_flds_mod, only: seq_flds_a2x_fields, seq_flds_x2a_fields !
  use seq_flds_mod, only: seq_flds_o2x_fields ! needed for MOAB init of ocean fields o2x to be able to transfer to coupler
  use seq_flds_mod, only: seq_flds_x2o_fields ! needed for MOAB init of ocean fields x2o to be able to transfer from coupler
  use seq_flds_mod, only: seq_flds_i2x_fields, seq_flds_x2i_fields ! needed for MOAB init of ice fields x2o on coupler side, to save them
  use seq_flds_mod, only: seq_flds_l2x_fields, seq_flds_x2l_fields !
  use seq_flds_mod, only: seq_flds_r2x_fields, seq_flds_x2r_fields, seq_flds_r2x_fluxes
  use seq_comm_mct, only: cplid, logunit
  use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs, seq_comm_iamin
  use seq_diag_mct
  use cplcomp_moab_helpers_mod, only: moab_register_app, moab_send_mesh, moab_receive_mesh, &
     moab_load_mesh, moab_free_sender_buffers, moab_define_global_id_tag, moab_define_double_tag

  use seq_comm_mct, only : mhid, mpoid, mbaxid, mboxid, mbofxid ! iMOAB app ids, for atm, ocean, ax mesh, ox mesh
  use seq_comm_mct, only : mhpgid         !    iMOAB app id for atm pgx grid, on atm pes
  use seq_comm_mct, only : atm_pg_active  ! flag if PG mesh instanced
  use seq_comm_mct, only : mlnid , mblxid !    iMOAB app id for land , on land pes and coupler pes
  use seq_comm_mct, only : mb_scm_land    !  logical used to identify land scm case; moab will migrate land then
  use seq_comm_mct, only : mb_dead_comps  !  logical to identify dead component configuration
  use seq_comm_mct, only : mphaid !            iMOAB app id for phys atm; comp atm is 5, phys 5+200
  use seq_comm_mct, only : MPSIID, mbixid  !  sea-ice on comp pes and on coupler pes
  use seq_comm_mct, only : mrofid, mbrxid  ! iMOAB id of moab rof app on comp pes and on coupler too
  use shr_mpi_mod,  only: shr_mpi_max
  ! use dimensions_mod, only : np     ! for atmosphere
  use iso_c_binding

  implicit none
  private  ! except
#include <mpif.h>
  save

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: cplcomp_moab_Init       ! called to migrate MOAB mesh from
                                    !   component pes to coupler pes
  public :: component_exch_moab
  public :: seq_mctext_avCreate
  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

   ! Shared routines for helper dispatch, invariants, and comm-graph setup.

   private :: copy_aream_from_area
   private :: moab_exchange_domain_tags
    private :: cplcomp_moab_init_atm
    private :: cplcomp_moab_init_ocn
    private :: cplcomp_moab_init_lnd
      private :: cplcomp_moab_init_ice
   private :: cplcomp_moab_init_rof
    private :: cplcomp_moab_resolve_comm_types
    private :: cplcomp_moab_compute_comm_graph
    private :: cplcomp_moab_atm_phys_cid
  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  integer,public :: seq_mctext_decomp

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(*),parameter :: subName = '(seq_mctext_mct)'
  real(r8),parameter :: c1 = 1.0_r8
   integer, parameter :: ATM_PHYS_ID_OFFSET = 200
   integer, parameter :: OCN_SECOND_COPY_ID_OFFSET = 1000

  !=======================================================================
contains
  !=======================================================================

  subroutine seq_mctext_avCreate(AVin,IDin,AVout,ID,lsize)

    !-----------------------------------------------------------------------
    ! Extend an AV to a larger set of pes or
    ! Initialize an AV on another set of pes
    !-----------------------------------------------------------------------

    implicit none
    type(mct_aVect), intent(INOUT):: AVin
    integer         ,intent(IN)   :: IDin ! ID associated with AVin
    type(mct_aVect), intent(INOUT):: AVout
    integer        , intent(IN)   :: ID   ! ID to initialize over
    integer        , intent(IN)   :: lsize

    ! Local variables

    character(len=*),parameter :: subname = "(seq_mctext_avCreate) "
    integer :: mpicom
    integer :: rank
    integer :: lsizei, lsizen
    integer :: srank,srankg
    integer :: ierr
    character(len=CXX) :: iList,rList

    call seq_comm_getinfo(ID,mpicom=mpicom,iam=rank)

    ! --- lsizen is the size of the newly initialized AV, zero is valid

    lsizei = -1
    if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
    lsizen = lsize

    ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
    ! --- set the pe and broadcast it to all other pes

    srank = -1
    srankg = -1
    if (lsizei > 0) srank = rank

    call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

    if (srankg < 0) then
       write(logunit,*) subname,' ERROR AVin not initialized '
       call shr_sys_abort()
    endif

    ! --- set the iList and rList from the broadcast pe (srankg) and
    ! --- broadcast the lists

    iList = " "
    rList = " "
    if (rank == srankg) then
       if (mct_aVect_nIAttr(AVin) /= 0) iList = mct_aVect_ExportIList2c(AVin)
       if (mct_aVect_nRattr(AVin) /= 0) rList = mct_aVect_ExportRList2c(AVin)
    endif

    call mpi_bcast(iList,len(iList),MPI_CHARACTER,srankg,mpicom,ierr)
    call mpi_bcast(rList,len(rList),MPI_CHARACTER,srankg,mpicom,ierr)

    ! --- now allocate the AV on all pes.  the AV should not exist before.
    ! --- If it does, mct should die.

    if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
       call mct_aVect_init(AVout,iList=iList,rList=rList,lsize=lsizen)
    elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
       call mct_aVect_init(AVout,iList=iList,lsize=lsizen)
    elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
       call mct_aVect_init(AVout,rList=rList,lsize=lsizen)
    endif

  end subroutine seq_mctext_avCreate

  !=======================================================================

subroutine  copy_aream_from_area(mbappid)

      ! maybe we will move this from here
      use iMOAB, only: iMOAB_GetDoubleTagStorage, iMOAB_SetDoubleTagStorage, iMOAB_GetMeshInfo

      integer, intent(in) :: mbappid
      character(CXX)           :: tagname
      integer                  nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
      real(r8),    allocatable    :: tagValues(:) ! used for setting aream tags for atm domain read case
      integer                     :: arrSize ! for the size of tagValues
      integer :: ierr, ent_type

      ! copy aream from area
      if (mbappid >= 0) then  ! coupler procs
         ierr  = iMOAB_GetMeshInfo ( mbappid, nvert, nvise, nbl, nsurf, nvisBC )
         ! Use mesh-intrinsic detection: if the mesh has cells, area
         ! lives on cells (domain-file meshes, PG2 ATM, dead comps).
         ! Point-cloud meshes (spectral ATM) have no cells, so use vertices.
         ! if (.not.atm_pg_active .and. .not. mb_dead_comps) then
         if (nvise(1) > 0) then
            arrSize  = nvise(1) ! cells
            ent_type = 1 ! cells
         else
            arrSize  = nvert(1) ! vertices (point cloud)
            ent_type = 0 ! vertices
         endif
         allocate(tagValues(arrSize))
         tagname = 'area'//C_NULL_CHAR
         ierr  = iMOAB_GetDoubleTagStorage( mbappid, tagname, arrsize , ent_type, tagValues )
         tagname = 'aream'//C_NULL_CHAR
         ierr  = iMOAB_SetDoubleTagStorage( mbappid, tagname, arrsize , ent_type, tagValues )
         deallocate(tagValues)
      endif

      return
  !=======================================================================
 end subroutine copy_aream_from_area

  integer function cplcomp_moab_atm_phys_cid(id_old)

      integer, intent(in) :: id_old

      cplcomp_moab_atm_phys_cid = ATM_PHYS_ID_OFFSET + id_old

  end function cplcomp_moab_atm_phys_cid

  subroutine moab_exchange_domain_tags(comp, comp_appid, cpl_appid, domain_fields, dom_context)
      type(component_type), intent(inout) :: comp
      integer,              intent(in)    :: comp_appid, cpl_appid
      character(len=*),     intent(in)    :: domain_fields, dom_context
      character(CXX) :: tagname
      tagname = trim(domain_fields)//C_NULL_CHAR
      call component_exch_moab(comp, comp_appid, cpl_appid, 'c2x', tagname, context_exch=dom_context)
      call copy_aream_from_area(cpl_appid)
  end subroutine moab_exchange_domain_tags

  subroutine cplcomp_moab_resolve_comm_types(src_has_cells, tgt_has_cells, typeA, typeB)

      logical, intent(in) :: src_has_cells, tgt_has_cells
      integer, intent(out) :: typeA, typeB

      if (src_has_cells) then
         typeA = 3
      else
         typeA = 2
      endif

      if (tgt_has_cells) then
         typeB = 3
      else
         typeB = 2
      endif

  end subroutine cplcomp_moab_resolve_comm_types

  subroutine cplcomp_moab_compute_comm_graph(src_appid, tgt_appid, mpicom_join, mpigrp_src, mpigrp_tgt, &
       src_has_cells, tgt_has_cells, src_context, tgt_context, subname, graph_name)

      use iMOAB, only: iMOAB_ComputeCommGraph

      integer, intent(in) :: src_appid, tgt_appid
      integer, intent(in) :: mpicom_join, mpigrp_src, mpigrp_tgt
      logical, intent(in) :: src_has_cells, tgt_has_cells
      integer, intent(in) :: src_context, tgt_context
      character(len=*), intent(in) :: subname, graph_name

      integer :: ierr, typeA, typeB

      if (mpicom_join == MPI_COMM_NULL) then
         write(logunit,*) subname,' invalid communicator for ', trim(graph_name)
         call shr_sys_abort(subname//' ERROR in computing comm graph for '//trim(graph_name))
      endif

      call cplcomp_moab_resolve_comm_types(src_has_cells, tgt_has_cells, typeA, typeB)

      ierr = iMOAB_ComputeCommGraph(src_appid, tgt_appid, mpicom_join, mpigrp_src, mpigrp_tgt, &
         typeA, typeB, src_context, tgt_context)
      if (ierr .ne. 0) then
         write(logunit,*) subname,' error in computing comm graph for ', trim(graph_name)
         call shr_sys_abort(subname//' ERROR in computing comm graph for '//trim(graph_name))
      endif

  end subroutine cplcomp_moab_compute_comm_graph

  subroutine cplcomp_moab_init_atm(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)

      use iMOAB, only: iMOAB_WriteMesh, iMOAB_GetMeshInfo
      use seq_infodata_mod, only: seq_infodata_type, seq_infodata_GetData

      type(seq_infodata_type), intent(in)  :: infodata
      type(component_type),    intent(inout) :: comp
      integer, intent(in) :: id_old, id_join
      integer, intent(in) :: mpicom_old, mpicom_new, mpicom_join
      logical, intent(in) :: dead_comps
      integer, intent(in) :: partMethod
      character(len=*), intent(in) :: subname

      integer :: mpigrp_cplid, mpigrp_old
      integer :: ierr
      character*200 :: appname, outfile, wopts, ropts, infile
      character(CL) :: atm_mesh
      integer :: ATM_PHYS_CID
      integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)

      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

      ! find atm mesh/domain file if it exists; it would be for data atm model (atm_prognostic false)
      call seq_infodata_GetData(infodata,atm_mesh = atm_mesh)

!!!!!!!! ON ATM COMPONENT
      if (mphaid >= 0) then  ! component atm procs
         ierr  = iMOAB_GetMeshInfo ( mphaid, nvert, nvise, nbl, nsurf, nvisBC )
         comp%mbApCCid = mphaid ! phys atm
         ! Auto-detect: dead FV ATM creates a full RLL mesh with cells,
         ! while spectral ATM sends only a point cloud (vertices).
         if (dead_comps) then
            comp%mbGridType = 1 ! cell mesh (dead FV ATM)
            comp%mblsize = nvise(1)
         else
            comp%mbGridType = 0 ! point cloud (spectral ATM)
            comp%mblsize = nvert(1)
         endif
      endif

!!!!!!!! ON ATM COMPONENT
      if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (atmosphere)
      !  send mesh to coupler
      !!!!  FULL ATM
         if ( trim(atm_mesh) == 'none' ) then ! full model
            if (atm_pg_active) then !  change : send the point cloud phys grid mesh, not coarse mesh,
                                    !     when atm pg active
               call moab_send_mesh(mhpgid, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
            else
               ! still use the mhid, original coarse mesh
               call moab_send_mesh(mphaid, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
            endif
         endif
      endif ! atmosphere pes
!!!!!!!!  ON ATM IN CPL
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         appname = "COUPLE_ATM"//C_NULL_CHAR
         ! migrated mesh gets another app id, moab atm to coupler (mbax)
         call moab_register_app(appname, mpicom_new, id_join, mbaxid, subname)
         !!!!  FULL ATM
         if ( trim(atm_mesh) == 'none' ) then ! full atm
            ! will receive either pg2 mesh, or point cloud mesh corresponding to GLL points
            ! (mphaid app) for spectral case
            ! this cannot be used for maps (either computed online or read)
            call moab_receive_mesh(mbaxid, mpicom_join, mpigrp_old, id_old, subname)
         !!!!  DATA ATM
         else
           ! we need to read the atm mesh on coupler, from domain file
            infile = trim(atm_mesh)//C_NULL_CHAR
            ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION;NO_CULLING'//C_NULL_CHAR
            if (seq_comm_iamroot(CPLID)) then
               write(logunit,'(A)') subname//' loading atm domain mesh from file '//trim(atm_mesh) &
                , ' with options ' // trim(ropts)
            endif
            call moab_load_mesh(mbaxid, infile, ropts, 0, subname)
            ! right now, turn atm_pg_active to true
            atm_pg_active = .true. ! FIXME TODO
            call moab_define_global_id_tag(mbaxid, subname)
         endif

      endif  ! on coupler pes
      !  iMOAB_FreeSenderBuffers needs to be called after receiving the mesh

!!!!!!!!  ATM COMPONENT
      if (mphaid .ge. 0) then  ! we are on component atm pes
!!!!! FULL ATM
         if ( trim(atm_mesh) == 'none' ) then  ! full atmosphere
            if (atm_pg_active) then! we send mesh from mhpgid app
               call moab_free_sender_buffers(mhpgid, id_join, subname)
            else
               ! we send mesh from point cloud data
               call moab_free_sender_buffers(mphaid, id_join, subname)
            endif
         endif
      endif  ! component atm pes

!!!!!  back to joint COMPONENT and CPL procs
      ! Graph between atm phys, mphaid, and atm dyn on coupler, mbaxid.
      ! Preserve the ATM physics offset invariant in one shared helper.
      ATM_PHYS_CID = cplcomp_moab_atm_phys_cid(id_old) ! 200 + 5 for atm, see line  969   ATM_PHYS = 200 + ATMID ! in
                                 !  components/cam/src/cpl/atm_comp_mct.F90
                                 !  components/data_comps/datm/src/atm_comp_mct.F90 ! line 177 !!

      ! this is not needed for migrating point cloud to point cloud !
      ! it is needed only after migrating pg2 mesh to cpupler
      call cplcomp_moab_compute_comm_graph(mphaid, mbaxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
          dead_comps, (atm_pg_active .or. dead_comps), ATM_PHYS_CID, id_join, subname, 'atm model')

      ! we can receive those tags only on coupler pes, when mbaxid exists
      ! we have to check that before we can define the tag
!!!!!! On ATM ON CPL
      if (mbaxid .ge. 0 ) then   !  coupler pes
         call moab_define_double_tag(mbaxid, trim(seq_flds_a2x_fields), subname)
         call moab_define_double_tag(mbaxid, trim(seq_flds_x2a_fields), subname)
         call moab_define_double_tag(mbaxid, trim(seq_flds_dom_fields)//":norm8wt", subname)
      endif ! coupler pes
      ! also, frac, area,  masks has to come from atm mphaid, not from domain file reader
      ! TODO:  this should be called on the joint procs, not coupler only.
      call moab_exchange_domain_tags(comp, mphaid, mbaxid, 'lat:lon:area:frac:mask', 'doma')

#ifdef MOABDEBUG
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         ! debug test
         if (atm_pg_active) then !
            outfile = 'recMeshAtmPG.h5m'//C_NULL_CHAR
         else
            outfile = 'recMeshAtm.h5m'//C_NULL_CHAR
         endif
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR
   !      write out the mesh file to disk
         ierr = iMOAB_WriteMesh(mbaxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing mesh '
            call shr_sys_abort(subname//' ERROR in writing mesh ')
         endif
      endif ! coupler pes
#endif
  end subroutine cplcomp_moab_init_atm

  subroutine cplcomp_moab_init_ocn(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)

      use iMOAB, only: iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo, iMOAB_ComputeCommGraph
      use seq_infodata_mod, only: seq_infodata_type, seq_infodata_GetData
      use shr_moab_mod, only: mbGetnCells, mbGetCellTagVals, mbSetCellTagVals

      type(seq_infodata_type), intent(in) :: infodata
      type(component_type),    intent(inout) :: comp
      integer, intent(in) :: id_old
      integer, intent(inout) :: id_join
      integer, intent(in) :: mpicom_old, mpicom_new, mpicom_join
      logical, intent(in) :: dead_comps
      integer, intent(in) :: partMethod
      character(len=*), intent(in) :: subname

      integer :: mpigrp_cplid, mpigrp_old
      integer :: ierr, context_id
      character*200 :: appname, outfile, wopts, ropts, infile
      character(CL) :: ocn_domain
      integer :: tagtype, numco, tagindex
      character(CXX) :: tagname
      integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
      real(r8), allocatable :: tagValues(:)
      integer :: arrsize, nloc

      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

      ! find ocean domain file if it exists; it would be for data ocean model (ocn_prognostic false)
      call seq_infodata_GetData(infodata,ocn_domain=ocn_domain)

!!!!!!  OCEAN COMPONENT
      if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (ocean)
#ifdef MOABDEBUG
         !   write out the mesh file to disk, in parallel
         !    we did it here because MOABDEBUG was not propagating with FFLAGS; we should move it
         !  now to component code, because MOABDEBUG can be propagated now with CPPDEFS
         outfile = 'wholeOcn.h5m'//C_NULL_CHAR
         wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
         ierr = iMOAB_WriteMesh(MPOID, outfile, wopts)
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing ocean mesh '
            call shr_sys_abort(subname//' ERROR in writing ocean mesh ')
         endif
#endif

         ierr  = iMOAB_GetMeshInfo ( mpoid, nvert, nvise, nbl, nsurf, nvisBC )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in getting mesh info on ocn '
            call shr_sys_abort(subname//' ERROR in getting mesh info on ocn  ')
         endif
         comp%mbApCCid = mpoid ! ocn comp app in moab
!!!!!  FULL OCN
         if ( trim(ocn_domain) == 'none' ) then
            !  send mesh to coupler
            call moab_send_mesh(mpoid, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
            comp%mbGridType = 1 ! cells
            comp%mblsize = nvise(1) ! cells
!!!!!  DATA OCN
         else
            comp%mbGridType = 0 ! vertices
            comp%mblsize = nvert(1) ! vertices
         endif
      endif
!!!!!!  OCN ON CPL
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         appname = "COUPLE_MPASO"//C_NULL_CHAR
         ! migrated mesh gets another app id, moab ocean to coupler (mbox)
         call moab_register_app(appname, mpicom_new, id_join, mboxid, subname)
 !!!!! FULL OCN
         if ( trim(ocn_domain) == 'none' ) then
            call moab_receive_mesh(mboxid, mpicom_join, mpigrp_old, id_old, subname)
 !!!!! DATA OCN
         else
           ! we need to read the ocean mesh on coupler, from domain file
            infile = trim(ocn_domain)//C_NULL_CHAR
            ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;NO_CULLING;REPARTITION'//C_NULL_CHAR
            if (seq_comm_iamroot(CPLID)) then
               write(logunit,'(A)') subname//' loading ocn domain mesh from file '//trim(infile) &
                , ' with options ' // trim(ropts)
            endif
            call moab_load_mesh(mboxid, infile, ropts, 0, subname)
            call moab_define_global_id_tag(mboxid, subname)
         endif  ! end of defining couplers copy of ocean mesh
 !!!! Still on OCN ON CPL
         tagname = trim(seq_flds_o2x_fields)//C_NULL_CHAR
         tagtype = 1  ! dense, double
         numco = 1 !  one value per cell
         ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags o2x on coupler'
            call shr_sys_abort(subname//' ERROR in defining tags o2x on coupler ')
         endif

         ! need also to define seq_flds_x2o_fields on coupler instance, and on ocean comp instance
         tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags x2o on coupler'
            call shr_sys_abort(subname//' ERROR in defining tags x2o on coupler ')
         endif

         !add the normalization tag
         tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on ocn on coupler '
            call shr_sys_abort(subname//' ERROR in defining tags ')
         endif

#ifdef MOABDEBUG
   !      debug test
         outfile = 'recMeshOcn.h5m'//C_NULL_CHAR
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
   !      write out the mesh file to disk
         ierr = iMOAB_WriteMesh(mboxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing ocean mesh coupler '
            call shr_sys_abort(subname//' ERROR in writing ocean mesh coupler ')
         endif
#endif
      endif
!!!!!!  OCEAN COMPONENT
      if (mpoid .ge. 0) then  ! we are on component ocn pes
         if ( trim(ocn_domain) == 'none' ) then
            context_id = id_join
            call moab_free_sender_buffers(mpoid, context_id, subname)
         endif
      endif
      ! in case of domain read, we need to compute the comm graph
!!!!!! on joint OCN and CPL procs
  !!!!! DATA OCN
      if ( trim(ocn_domain) /= 'none' ) then
         ! we are now on joint pes, compute comm graph between data ocn and coupler model ocn
         call cplcomp_moab_compute_comm_graph(mpoid, mboxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
            dead_comps, .true., id_old, id_join, subname, 'data ocn model')
         ! also, frac, area,  masks has to come from ocean mpoid, not from domain file reader
         ! this is hard to digest :(
         tagname = 'lat:lon:area:frac:mask'//C_NULL_CHAR
         call component_exch_moab(comp, mpoid, mboxid, 'c2x', tagname, context_exch='domo')
         if (mboxid > 0) then ! on coupler pes only
            call copy_aream_from_area(mboxid)
         endif
      endif

!!!!!!!!! OCEAN 2nd COPY
      ! start copy
      ! Do another ocean copy on the coupler so So_fswpen does not collide between
      ! xao states and o2x states. Preserve the explicit second-copy context offset.
      id_join = id_join + OCN_SECOND_COPY_ID_OFFSET

!!!!!!  OCEAN COMPONENT
      if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (ocean)
 !!!! FULL OCEAN
         if ( trim(ocn_domain) == 'none' ) then
            !  send mesh to coupler, the second time! a copy would be cheaper
            call moab_send_mesh(mpoid, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
         endif
      endif

!!!!!!  ON 2nd OCN ON CPL
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         appname = "COUPLE_MPASOF"//C_NULL_CHAR
         ! migrated mesh gets another app id, moab ocean to coupler (mbox)
         call moab_register_app(appname, mpicom_new, id_join, mbofxid, subname)
 !!!!! FULL OCN
         if ( trim(ocn_domain) == 'none' ) then
            call moab_receive_mesh(mbofxid, mpicom_join, mpigrp_old, id_old, subname)
 !!!!! DATA OCN
         else
             ! we need to read the ocean mesh on coupler, from domain file
            infile = trim(ocn_domain)//C_NULL_CHAR
            ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;NO_CULLING;REPARTITION'//C_NULL_CHAR
            if (seq_comm_iamroot(CPLID)) then
               write(logunit,'(A)') subname//' load ocn domain mesh from file for second ocn instance '//trim(ocn_domain) &
               , ' with options '//trim(ropts)
            endif
            call moab_load_mesh(mbofxid, infile, ropts, 0, subname)
            call moab_define_global_id_tag(mbofxid, subname)
         endif

  !! ON 2nd OCEAN ON COUPLER
         tagtype = 1  ! dense, real
         numco = 1
         tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbofxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on ocn on coupler '
            call shr_sys_abort(subname//' ERROR in defining tags ')
         endif

         ! copy domain data to mbofxid
         tagname = 'lat:lon:area:frac:mask'//C_NULL_CHAR
         nloc = mbGetnCells(mbofxid)
         arrsize=nloc*5
         allocate(tagValues(arrsize))
         call mbGetCellTagVals(mboxid,tagname,tagValues,arrsize)
         call mbSetCellTagVals(mbofxid,tagname,tagValues,arrsize)
         deallocate(tagValues)

      endif

!!!!!!  ON OCN COMPONENT
      if (mpoid .ge. 0) then  ! we are on component ocn pes again, release buffers
          if ( trim(ocn_domain) == 'none' ) then
            context_id = id_join
                   call moab_free_sender_buffers(mpoid, context_id, subname)
          endif
      endif
      ! end copy
#ifdef MOABDEBUG
   if (mbofxid >= 0) then
      outfile = 'recMeshOcnF.h5m'//C_NULL_CHAR
      wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
   !  write out the mesh file to disk
      ierr = iMOAB_WriteMesh(mbofxid, trim(outfile), trim(wopts))
      if (ierr .ne. 0) then
          write(logunit,*) subname,' error in writing ocean mesh coupler '
          call shr_sys_abort(subname//' ERROR in writing ocean mesh coupler ')
      endif
   endif
#endif

  end subroutine cplcomp_moab_init_ocn

  subroutine cplcomp_moab_init_lnd(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)

      use iMOAB, only: iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorage, iMOAB_ComputeCommGraph
      use seq_infodata_mod, only: seq_infodata_type, seq_infodata_GetData

      type(seq_infodata_type), intent(in) :: infodata
      type(component_type),    intent(inout) :: comp
      integer, intent(in) :: id_old, id_join
      integer, intent(in) :: mpicom_old, mpicom_new, mpicom_join
      logical, intent(in) :: dead_comps
      integer, intent(in) :: partMethod
      character(len=*), intent(in) :: subname

      integer :: mpigrp_cplid, mpigrp_old
      integer :: ierr, context_id
      character*200 :: appname, outfile, wopts, ropts
      character(CL) :: lnd_domain
      integer :: tagtype, numco, tagindex, nghlay
      integer :: ent_type
      character(CXX) :: tagname, newlist
      integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
      logical :: rof_present, lnd_prognostic, single_column, scm_multcols
      real(r8), allocatable :: tagValues(:)
      integer :: arrsize, nfields
      type(mct_list) :: temp_list

      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes
      call seq_infodata_GetData(infodata,rof_present=rof_present, lnd_prognostic=lnd_prognostic)
      call seq_infodata_GetData(infodata,single_column=single_column, &
	 scm_multcols=scm_multcols)
      if (single_column .or. scm_multcols) then
         ! SCM land uses migrated mesh; do not read from domain file.
         mb_scm_land = .true.
      endif
      call seq_infodata_GetData(infodata,lnd_domain=lnd_domain)

      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         appname = "COUPLE_LAND"//C_NULL_CHAR
         ! migrated mesh gets another app id, moab land to coupler (mblx)
         call moab_register_app(appname, mpicom_new, id_join, mblxid, subname)
      endif

      if (mb_scm_land .or. dead_comps) then
          !  change : send the point cloud land mesh (1 point usually)
          !     when mb_scm_land
         if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (land)
            !  send mesh to coupler then
            call moab_send_mesh(mlnid, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
         endif
         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            call moab_receive_mesh(mblxid, mpicom_join, mpigrp_old, id_old, subname)
         endif
         if (MPI_COMM_NULL /= mpicom_old) then  ! we are on component lnd pes again, release buffers
            context_id = id_join
            call moab_free_sender_buffers(mlnid, context_id, subname)
        endif
      else
         ! Coupler side only: read land mesh from domain file
         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
               ! do not cull in case of data land, like all other data models
               ! for regular land model, cull, because the lnd component culls too
               if (lnd_prognostic) then
                  ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION'//C_NULL_CHAR
               else
                  ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION;NO_CULLING'//C_NULL_CHAR
               endif
               outfile = trim(lnd_domain)//C_NULL_CHAR
               nghlay = 0 ! no ghost layers
               if (seq_comm_iamroot(CPLID) ) then
                  write(logunit, *) "loading land domain file from file: ", trim(lnd_domain), &
                  " with options: ", trim(ropts)
               endif
               call moab_load_mesh(mblxid, outfile, ropts, nghlay, subname)
         endif ! end of coupler pes
      endif ! end of mb_scm_land check

      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         call moab_define_global_id_tag(mblxid, subname)

!  need to define tags on land too
         tagname = trim(seq_flds_l2x_fields)//C_NULL_CHAR
         tagtype = 1  ! dense, double
         numco = 1 !  one value per cell
         ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco, tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags l2x on coupler land'
            call shr_sys_abort(subname//' ERROR in defining tags l2x on coupler ')
         endif
         ! need also to define seq_flds_x2l_fields on coupler instance, and on land comp instance
         tagname = trim(seq_flds_x2l_fields)//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco, tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags x2l on coupler land'
            call shr_sys_abort(subname//' ERROR in defining tags x2l on coupler land')
         endif

         if (.not.rof_present) then  ! need to zero out some Flrr fields
            call shr_string_listIntersect(seq_flds_x2l_fields,seq_flds_r2x_fluxes,newlist)
            call mct_list_init(temp_list, newlist)
            nfields=mct_list_nitem (temp_list)
            if (nfields > 0) then
              ierr  = iMOAB_GetMeshInfo ( mblxid, nvert, nvise, nbl, nsurf, nvisBC )
              if (mb_scm_land) then
                 arrsize = nvert(1)*nfields
                 ent_type = 0 ! cell
              else
                 arrsize = nvise(1)*nfields
                 ent_type = 1 ! cell
              endif
              allocate(tagValues(arrsize))
              tagname = trim(newlist)//C_NULL_CHAR
              tagValues = 0.0_r8
              ierr = iMOAB_SetDoubleTagStorage ( mblxid, tagname, arrsize, ent_type, tagValues)
              if (ierr .ne. 0) then
                 write(logunit,*) subname,' error in zeroing Flrr tags on land', ierr
                 call shr_sys_abort(subname//' ERROR in zeroing Flrr tags land')
              endif
              deallocate(tagValues)
            endif
         endif

         !add the normalization tag
         tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on lnd on coupler '
            call shr_sys_abort(subname//' ERROR in defining tags ')
         endif

      endif ! end of coupler pes

      if (mlnid >= 0) then
         ierr  = iMOAB_GetMeshInfo ( mlnid, nvert, nvise, nbl, nsurf, nvisBC )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in getting mesh info for lnd on coupler '
            call shr_sys_abort(subname//' ERROR in getting mesh info for lnd on coupler ')
         endif
         comp%mbApCCid = mlnid ! land
         if (dead_comps) then
            comp%mbGridType = 1 ! dead comps create full mesh
            comp%mblsize = nvise(1) ! cells
         else
            comp%mbGridType = 0 ! 0 or 1, pc or cells
            comp%mblsize = nvert(1) ! vertices
         endif
      endif
      if ( .not. mb_scm_land ) then
         ! compute comm graph for tag exchange between land component and coupler
         call cplcomp_moab_compute_comm_graph(mlnid, mblxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
            dead_comps, .true., id_old, id_join, subname, 'lnd model')
      endif
      tagname = 'lat:lon:area:frac:mask'//C_NULL_CHAR
      call component_exch_moab(comp, mlnid, mblxid, 'c2x', tagname, context_exch='doml')
      if (mblxid > 0) then ! on coupler pes only
         call copy_aream_from_area(mblxid)
      endif

#ifdef MOABDEBUG
         outfile = 'recMeshLand.h5m'//C_NULL_CHAR
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
   !       write out the mesh file to disk
         ierr = iMOAB_WriteMesh(mblxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing land coupler mesh'
            call shr_sys_abort(subname//' ERROR in writing land coupler mesh')
         endif
#endif

  end subroutine cplcomp_moab_init_lnd

  subroutine cplcomp_moab_init_ice(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)

      use iMOAB, only: iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo, iMOAB_ComputeCommGraph
      use seq_infodata_mod, only: seq_infodata_type, seq_infodata_GetData

      type(seq_infodata_type), intent(in) :: infodata
      type(component_type),    intent(inout) :: comp
      integer, intent(in) :: id_old, id_join
      integer, intent(in) :: mpicom_old, mpicom_new, mpicom_join
      logical, intent(in) :: dead_comps
      integer, intent(in) :: partMethod
      character(len=*), intent(in) :: subname

      integer :: mpigrp_cplid, mpigrp_old
      integer :: ierr, context_id
      character*200 :: appname, outfile, wopts, ropts, infile
      character(CL) :: ice_domain
      integer :: tagtype, numco, tagindex
      character(CXX) :: tagname
      integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)

      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes
      ! find ice domain file if it exists; it would be for data ice model (ice_prognostic false)
      call seq_infodata_GetData(infodata,ice_domain=ice_domain)
      if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component p
#ifdef MOABDEBUG
         outfile = 'wholeSeaIce.h5m'//C_NULL_CHAR
         wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
         ierr = iMOAB_WriteMesh(MPSIID, outfile, wopts)
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing sea-ice'
            call shr_sys_abort(subname//' ERROR in writing sea-ice')
         endif
#endif
! start copy from ocean code
         if (MPSIID >= 0) then
            ierr  = iMOAB_GetMeshInfo ( MPSIID, nvert, nvise, nbl, nsurf, nvisBC )
            comp%mbApCCid = MPSIID ! ice imoab app id
         endif
         if ( trim(ice_domain) == 'none' ) then ! regular ice model
            if (dead_comps) then
               comp%mbGridType = 1 ! dead comps create full mesh
               comp%mblsize = nvise(1) ! cells
            else
               comp%mbGridType = 1 ! 0 or 1, pc or cells
               comp%mblsize = nvise(1) ! cells
            endif
            !  send sea ice mesh to coupler
            call moab_send_mesh(MPSIID, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
         else
            ! we could be using cice model
            comp%mbGridType = 0 ! 0 or 1, pc or cells
            comp%mblsize = nvert(1) ! vertices
         endif
      endif
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         appname = "COUPLE_MPASSI"//C_NULL_CHAR
         ! migrated mesh gets another app id, moab moab sea ice to coupler (mbix)
         call moab_register_app(appname, mpicom_new, id_join, mbixid, subname)
         if ( trim(ice_domain) == 'none' ) then ! regular ice model
            call moab_receive_mesh(mbixid, mpicom_join, mpigrp_old, id_old, subname)
         else
            ! we need to read the mesh ice (domain file)
            ! we could be using cice model or data sea ice; in both cases ice_domain should be non-empty
            ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;NO_CULLING;REPARTITION'//C_NULL_CHAR
            infile = trim(ice_domain)//C_NULL_CHAR
            if (seq_comm_iamroot(CPLID)) then
               write(logunit,'(A)') subname//' loading ice domain mesh from file '//infile &
                 , ' with options '//trim(ropts)
            endif
            call moab_load_mesh(mbixid, infile, ropts, 0, subname)
            call moab_define_global_id_tag(mbixid, subname)
         endif ! end data ice

         if (MPSIID .ge. 0) then  ! we are on component sea ice pes
            if ( trim(ice_domain) == 'none' ) then
               context_id = id_join
               call moab_free_sender_buffers(MPSIID, context_id, subname)
            endif
         endif

         tagtype = 1  ! dense, double
         numco = 1 !  one value per cell / entity
         tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
         if ( ierr == 1 ) then
            call shr_sys_abort( subname//' ERROR: cannot define tags for ice on coupler' )
         end if
         tagname = trim(seq_flds_x2i_fields)//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
         if ( ierr == 1 ) then
            call shr_sys_abort( subname//' ERROR: cannot define tags for ice on coupler' )
         end if

         !add the normalization tag
         tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on ice on coupler '
            call shr_sys_abort(subname//' ERROR in defining tags ')
         endif

         ! add data that is interpolated to sea ice
         tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
         tagtype = 1 ! dense
         numco = 1 !
         ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags for seq_flds_a2x_fields on ice cpl'
            call shr_sys_abort(subname//' ERROR in coin defining tags for seq_flds_a2x_fields on ice cpl')
         endif

         ! add data that is interpolated to sea ice
         tagname = trim(seq_flds_r2x_fields)//C_NULL_CHAR
         tagtype = 1 ! dense
         numco = 1 !
         ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags for seq_flds_r2x_fields on ice cpl'
            call shr_sys_abort(subname//' ERROR in coin defining tags for seq_flds_a2x_fields on ice cpl')
         endif

      endif

     ! in case of ice domain read, we need to compute the comm graph
     if ( trim(ice_domain) /= 'none' ) then
         ! we are now on joint pes, compute comm graph between data ice and coupler model ice
         call cplcomp_moab_compute_comm_graph(MPSIID, mbixid, mpicom_join, mpigrp_old, mpigrp_cplid, &
            dead_comps, .true., id_old, id_join, subname, 'data ice model')
         ! also, frac, area,  masks has to come from ice MPSIID , not from domain file reader
         ! this is hard to digest :(
         tagname = 'lat:lon:area:frac:mask'//C_NULL_CHAR
         call component_exch_moab(comp, MPSIID, mbixid, 'c2x', tagname, context_exch='domi')

         if (mbixid > 0) then ! on coupler pes only
            call copy_aream_from_area(mbixid)
         endif
      endif
#ifdef MOABDEBUG
!      debug test
      outfile = 'recMeshSeaIce.h5m'//C_NULL_CHAR
      wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
!      write out the mesh file to disk
      ierr = iMOAB_WriteMesh(mbixid, trim(outfile), trim(wopts))
      if (ierr .ne. 0) then
          write(logunit,*) subname,' error in writing sea ice mesh on coupler '
          call shr_sys_abort(subname//' ERROR in writing sea ice mesh on coupler ')
      endif
#endif

  end subroutine cplcomp_moab_init_ice

  subroutine cplcomp_moab_init_rof(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)

      use iMOAB, only: iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo, iMOAB_ComputeCommGraph
      use seq_infodata_mod, only: seq_infodata_type, seq_infodata_GetData

      type(seq_infodata_type), intent(in) :: infodata
      type(component_type),    intent(inout) :: comp
      integer, intent(in) :: id_old, id_join
      integer, intent(in) :: mpicom_old, mpicom_new, mpicom_join
      logical, intent(in) :: dead_comps
      integer, intent(in) :: partMethod
      character(len=*), intent(in) :: subname

      integer :: mpigrp_cplid, mpigrp_old
      integer :: ierr, context_id
      character*200 :: appname, outfile, wopts, ropts
      character(CL) :: rtm_mesh, rof_domain
      integer :: tagtype, numco, tagindex, nghlay
      character(CXX) :: tagname
      integer :: nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)

      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes
      call seq_infodata_GetData(infodata,rof_domain=rof_domain,rof_mesh=rtm_mesh)

      ! Component side: send mesh to coupler (BEFORE coupler registers and receives)
      if (mrofid >= 0) then  ! component pes
         ierr  = iMOAB_GetMeshInfo ( mrofid, nvert, nvise, nbl, nsurf, nvisBC )
         comp%mbApCCid = mrofid !
         if (dead_comps) then
            comp%mbGridType = 1 ! dead comps create full mesh
            comp%mblsize = nvise(1) ! cells
         else
            comp%mbGridType = 0 ! 0 or 1, pc or cells
            comp%mblsize = nvert(1) ! vertices
         endif
      endif

      if (dead_comps) then ! full river model
         if (MPI_COMM_NULL /= mpicom_old .and. mrofid >= 0) then ! component pes, send mesh to coupler
            call moab_send_mesh(mrofid, mpicom_join, mpigrp_cplid, id_join, partMethod, subname)
         endif ! component pes
      endif

      ! Coupler side: register application and receive/load mesh
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
         appname = "COUPLE_MROF"//C_NULL_CHAR
         call moab_register_app(appname, mpicom_new, id_join, mbrxid, subname)

         if (dead_comps) then
            ! migrated mesh gets another app id, moab rof to coupler (mbrx)
            ! Receive mesh from river component
            call moab_receive_mesh(mbrxid, mpicom_join, mpigrp_old, id_old, subname)
         else
            ! we will read the mesh from domain file
            ! first check if we already have elements in mbrxid
            ierr  = iMOAB_GetMeshInfo ( mbrxid, nvert, nvise, nbl, nsurf, nvisBC )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in getting mesh info on ROF coupler '
               call shr_sys_abort(subname//' ERROR in getting mesh info on ROF coupler  ')
            endif
            if (nvert(3) .eq. 0) then
               ! load mesh from scrip file passed from river model, if domain file is not available
               call seq_infodata_GetData(infodata,rof_mesh=rtm_mesh,rof_domain=rof_domain)
               if ( trim(rof_domain) == 'none' ) then
                  outfile = trim(rtm_mesh)//C_NULL_CHAR
                  ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN'//C_NULL_CHAR
               else
                  outfile = trim(rof_domain)//C_NULL_CHAR
                  ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION'//C_NULL_CHAR
               endif
               nghlay = 0 ! no ghost layers
               if (seq_comm_iamroot(CPLID)) then
                  write(logunit,'(A)') subname//' loading rof from file '//trim(outfile) &
                     , ' with options ', trim(ropts)
               endif
               call moab_load_mesh(mbrxid, outfile, ropts, nghlay, subname)
            endif ! end of reading rof mesh from file

         end if ! end of dead_comps check

         call moab_define_global_id_tag(mbrxid, subname)

         tagtype = 1  ! dense, double
         numco = 1 !  one value per cell / entity
         tagname = trim(seq_flds_r2x_fields)//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco, tagindex )
         if ( ierr == 1 ) then
            call shr_sys_abort( subname//' ERROR: cannot define tags for rof on coupler' )
         end if
         tagname = trim(seq_flds_x2r_fields)//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco, tagindex )
         if ( ierr == 1 ) then
            call shr_sys_abort( subname//' ERROR: cannot define tags for rof on coupler' )
         end if

         !add the normalization tag
         tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
         ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on rof on coupler '
            call shr_sys_abort(subname//' ERROR in defining tags ')
         endif

      endif  ! coupler pes

      ! Free buffers on component side after send completes
      if (mrofid >= 0 .and. dead_comps) then
         context_id = id_join
         call moab_free_sender_buffers(mrofid, context_id, subname)
      endif

      ! we are now on joint pes, compute comm graph between rof and coupler model
      ! typeA=3 for dead comps (full mesh), typeA=2 for regular rof (point cloud)
      ! typeB=3 always: coupler always has full mesh (loaded from file or received from dead comp)
      call cplcomp_moab_compute_comm_graph(mrofid, mbrxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
         dead_comps, .true., id_old, id_join, subname, 'rof model')

      tagname = 'area:lon:lat:frac:mask'//C_NULL_CHAR
      call component_exch_moab(comp, mrofid, mbrxid, 'c2x', tagname, context_exch='domr')
      ! copy aream from area in all cases
      ! initialize aream from area; it may have different values in the end, or reset again
      if (mbrxid > 0) then ! on coupler pes only
         call copy_aream_from_area(mbrxid)
      endif
#ifdef MOABDEBUG
      if (mbrxid >= 0) then
         outfile = 'recMeshRof.h5m'//C_NULL_CHAR
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
!         write out the mesh file to disk
         ierr = iMOAB_WriteMesh(mbrxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
           write(logunit,*) subname,' error in writing rof mesh on coupler '
            call shr_sys_abort(subname//' ERROR in writing rof mesh on coupler ')
         endif
      endif
#endif

  end subroutine cplcomp_moab_init_rof

   subroutine cplcomp_moab_Init(infodata,comp)

      ! This routine initializes iMOAB applications for each component,
      !  dispatching to component-specific helpers for mesh migration,
      !  tag setup, and comm-graph construction on coupler/component PEs.

      !-----------------------------------------------------
      !
      use iMOAB, only: iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo, &
      iMOAB_SetDoubleTagStorage, iMOAB_ComputeCommGraph
      !
      use seq_infodata_mod
      use shr_moab_mod
      !
      type(seq_infodata_type) ,  intent(in) :: infodata
      type(component_type), intent(inout) :: comp

      !
      ! Local Variables
      !
      integer                  :: mpicom_cplid
      integer                  :: mpicom_old
      integer                  :: mpicom_new
      integer                  :: mpicom_join
      integer                  :: ID_old
      integer                  :: ID_new
      integer                  :: ID_join

      character(len=*),parameter :: subname = "(cplcomp_moab_Init) "

      integer                  :: mpigrp_cplid ! coupler pes
      integer                  :: mpigrp_old   !  component group pes
      integer                  :: ierr, context_id
      character*200            :: appname, outfile, wopts, ropts, infile
      character(CL)            :: rtm_mesh, rof_domain
      character(CL)            :: lnd_domain
      character(CL)            :: ocn_domain
      character(CL)            :: ice_domain   ! used for data ice only?
      character(CL)            :: atm_mesh
      integer                  :: maxMH, maxMPO, maxMLID, maxMSID, maxMRID ! max pids for moab apps atm, ocn, lnd, sea-ice, rof
      integer                  :: tagtype, numco,  tagindex, partMethod, nghlay
      integer                  :: rank, ent_type
      integer                  :: typeA, typeB, ATM_PHYS_CID ! used to compute par graph between atm phys
                                                            ! and atm spectral on coupler
      character(CXX)           :: tagname
      character(CXX)           :: newlist
      integer                  nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
      logical                  :: rof_present, lnd_prognostic, single_column, scm_multcols
      real(r8),    allocatable    :: tagValues(:) ! used for setting aream tags for atm domain read case
      integer                     :: arrsize ! for the size of tagValues
      type(mct_list)             :: temp_list
      integer                    :: nfields,nloc
      logical                    :: dead_comps
      ! real(R8), allocatable, target :: values(:)


   !-----------------------------------------------------

      call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID, iam=rank)

      id_new  = CPLID
      id_old  = comp%compid
      id_join = comp%cplcompid

      mpicom_new  = mpicom_CPLID
      mpicom_old  = comp%mpicom_compid
      mpicom_join = comp%mpicom_cplcompid

      ! partMethod = 0 ! trivial partitioning
      partMethod = 2 ! it is better to use RCB for atmosphere and ocean (needs MOAB_HAVE_ZOLTAN)
      context_id = -1 ! original sends/receives, so the context is -1
                     ! needed only to free send buffers

      call seq_comm_getinfo(ID_old ,mpicom=mpicom_old)
      call seq_comm_getinfo(ID_new ,mpicom=mpicom_new)
      call seq_comm_getinfo(ID_join,mpicom=mpicom_join)

      call shr_mpi_max(mphaid, maxMH, mpicom_join, all=.true.) ! if on atm / cpl joint, maxMH /= -1
      call shr_mpi_max(mpoid, maxMPO, mpicom_join, all=.true.)
      call shr_mpi_max(mlnid, maxMLID, mpicom_join, all=.true.)
      call shr_mpi_max(MPSIID, maxMSID, mpicom_join, all=.true.)
      call shr_mpi_max(mrofid, maxMRID, mpicom_join, all=.true.)

      if (seq_comm_iamroot(CPLID) ) then
         write(logunit, *) "MOAB coupling for ", comp%oneletterid,' ', comp%ntype
      endif

      call seq_infodata_GetData(infodata,dead_comps=dead_comps)
      mb_dead_comps = dead_comps

!!!!!!!!!!!!!!!! ATMOSPHERE
      if ( comp%oneletterid == 'a' .and. maxMH /= -1) then
         call cplcomp_moab_init_atm(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)
      endif  ! comp%oneletterid == 'a'

!!!!!!!!!!!!!!!! OCEAN
      ! ocean
      if (comp%oneletterid == 'o'  .and. maxMPO /= -1) then
         call cplcomp_moab_init_ocn(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)
      endif  ! end ocean

!!!!!!!!!!!!!!!! LAND
      if (comp%oneletterid == 'l'  .and. maxMLID /= -1) then
         call cplcomp_moab_init_lnd(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)
      endif  ! End of land model

!!!!!!!!!!!!!!!! SEA ICE
      if (comp%oneletterid == 'i'  .and. maxMSID /= -1) then
         call cplcomp_moab_init_ice(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)
      endif

!!!!!!!!!!!!!!!! RIVER
      if (comp%oneletterid == 'r'  .and. maxMRID /= -1) then
         call cplcomp_moab_init_rof(infodata, comp, id_old, id_join, mpicom_old, mpicom_new, mpicom_join, dead_comps, partMethod, subname)
      endif ! end for rof coupler set up

   end subroutine cplcomp_moab_Init


  !===============================================================================
  ! component_exch_moab
  !
  ! PURPOSE:
  !   Exchange field data (tags) between component and coupler meshes using iMOAB.
  !   This routine handles bidirectional data transfer in the first hop of a 2-hop
  !   exchange pattern for MOAB-based coupling.
  !
  ! DESCRIPTION:
  !   This subroutine coordinates MPI-based field data exchange between MOAB mesh
  !   applications on component PEs and coupler PEs. It uses iMOAB's send/receive
  !   infrastructure to transfer tagged data, manages communication buffers, and
  !   provides optional timing and debugging capabilities.
  !
  ! FLOW:
  !   1. Execute optional barrier timing
  !   2. Start exchange and map timers
  !   3. Determine source/target IDs based on direction (c2x or x2c)
  !   4. Send field tags from source mesh application
  !   5. Receive field tags at target mesh application
  !   6. Free sender communication buffers
  !   7. Stop map timer
  !   8. Exchange infodata if provided
  !   9. Stop exchange timer
  !   10. Optional debug output (if MOABDEBUG defined)
  !
  ! ARGUMENTS:
  !   comp              - component type containing MPI communicator info
  !   mbAPPid1          - iMOAB application ID for first mesh (sender or receiver)
  !   mbAppid2          - iMOAB application ID for second mesh (receiver or sender)
  !   direction         - data flow direction: 'c2x' (component->coupler) or
  !                       'x2c' (coupler->component)
  !   fields            - colon-separated list of field names to exchange
  !   context_exch      - optional context string for debugging output
  !   infodata          - optional metadata exchange object
  !   infodata_string   - optional string for infodata exchange
  !   mpicom_barrier    - optional MPI communicator for barriers
  !   run_barriers      - optional flag to enable/disable barriers
  !   timer_barrier     - optional timer name for barrier timing
  !   timer_comp_exch   - optional timer name for component exchange
  !   timer_map_exch    - optional timer name for mapping exchange
  !   timer_infodata_exch - optional timer name for infodata exchange
  !
  ! NOTES:
  !   - For atmosphere component, component-side ID is adjusted by +200 to handle
  !     point cloud representation
  !   - Sender buffers are freed after data transfer to conserve memory
  !   - Debug output writes mesh files when MOABDEBUG is defined
  !
  !===============================================================================
  subroutine component_exch_moab(comp, mbAPPid1, mbAppid2, direction, fields, context_exch, &
       infodata, infodata_string, mpicom_barrier, run_barriers, &
       timer_barrier, timer_comp_exch, timer_map_exch, timer_infodata_exch)

   use iMOAB ,  only: iMOAB_SendElementTag, iMOAB_ReceiveElementTag, iMOAB_WriteMesh, iMOAB_FreeSenderBuffers
   use seq_comm_mct, only :  num_moab_exports ! for debugging
   use ISO_C_BINDING, only : C_NULL_CHAR
   use shr_kind_mod      , only :  CXX => shr_kind_CXX
   use seq_infodata_mod, only: seq_infodata_exchange, seq_infodata_type
   use t_drv_timers_mod
   !---------------------------------------------------------------

    type(component_type)     , intent(in)           :: comp
    ! direction 'c2x' is from component to coupler; 'x2c' is from coupler to component
    integer,                   intent(in)           :: mbAPPid1, mbAppid2
    character(len=*)         , intent(in)           :: direction
    character(CXX)           , intent(in)           :: fields
    character(len=*)        ,  intent(in), optional :: context_exch
    type(seq_infodata_type) , intent(inout), optional :: infodata
    character(len=*)        , intent(in), optional :: infodata_string
    integer                 , intent(in), optional :: mpicom_barrier
    logical                 , intent(in), optional :: run_barriers
    character(len=*)        , intent(in), optional :: timer_barrier
    character(len=*)        , intent(in), optional :: timer_comp_exch
    character(len=*)        , intent(in), optional :: timer_map_exch
    character(len=*)        , intent(in), optional :: timer_infodata_exch

    character(*), parameter :: subname = '(component_exch_moab)'
    integer :: id_join, source_id, target_id, ierr
    integer :: mpicom_join
    character(CXX)              :: tagname
    character*100 outfile, wopts, lnum

    !---------------------------------------------------------------------------
    ! Optional barrier timing for performance analysis
    !---------------------------------------------------------------------------
    if (present(timer_barrier)) then
       if (present(run_barriers)) then
          if (run_barriers) then
             call t_drvstartf (trim(timer_barrier))
             call mpi_barrier(comp%mpicom_cplallcompid, ierr)
             call t_drvstopf (trim(timer_barrier))
          endif
       endif
    end if

    !---------------------------------------------------------------------------
    ! Start performance timers for component exchange and mapping
    !---------------------------------------------------------------------------
    if (present(timer_comp_exch)) then
      if (present(mpicom_barrier)) then
         call t_drvstartf (trim(timer_comp_exch), cplcom=.true., barrier=mpicom_barrier)
      end if
    end if

    if (comp%iamin_cplcompid) then
       if (present(timer_map_exch)) then
          call t_drvstartf (trim(timer_map_exch), barrier=comp%mpicom_cplcompid)
       end if

       !---------------------------------------------------------------------------
       ! Get joint communicator spanning both component and coupler PEs
       !---------------------------------------------------------------------------
       id_join = comp%cplcompid
       call seq_comm_getinfo(ID_join,mpicom=mpicom_join)

       ! Prepare tag name with C null terminator for iMOAB interface
       tagName = trim(fields)//C_NULL_CHAR

       !---------------------------------------------------------------------------
       ! Determine source and target IDs based on data flow direction
       ! c2x: component to coupler, x2c: coupler to component
       !---------------------------------------------------------------------------
       if (direction .eq. 'c2x') then
          source_id = comp%compid
          target_id = comp%cplcompid
       else if (direction .eq. 'x2c') then
          source_id = comp%cplcompid
          target_id = comp%compid
       else
          call shr_sys_abort(subname//' invalid direction in component_exch_moab: '// &
     &                         trim(direction))
       endif

       !---------------------------------------------------------------------------
       ! Special handling for atmosphere: add 200 to component-side ID
       ! This offset accounts for the point cloud representation used in
       ! atmosphere physics (vs spectral dynamics). The +200 convention matches
       ! the ATM_PHYS_CID offset used elsewhere in the coupler code.
       !---------------------------------------------------------------------------
       if (comp%oneletterid == 'a' .and. direction .eq. 'c2x' ) then
          source_id = source_id + ATM_PHYS_ID_OFFSET
       endif
       if (comp%oneletterid == 'a' .and. direction .eq. 'x2c' ) then
          target_id = target_id + ATM_PHYS_ID_OFFSET
       endif

       !---------------------------------------------------------------------------
       ! Send field tags from source mesh application
       ! Only PEs with valid mbAPPid1 participate in sending
       !---------------------------------------------------------------------------
       if (mbAPPid1 .ge. 0) then !  we are on the sending pes
          ierr = iMOAB_SendElementTag(mbAPPid1, tagName, mpicom_join, target_id)
          if (ierr .ne. 0) then
             call shr_sys_abort(subname//' cannot send element tag: '//trim(tagName))
          endif
       endif

       !---------------------------------------------------------------------------
       ! Receive field tags at target mesh application
       ! Only PEs with valid mbAPPid2 participate in receiving
       !---------------------------------------------------------------------------
       if ( mbAPPid2 .ge. 0 ) then !  we are on receiving end
          ierr = iMOAB_ReceiveElementTag(mbAPPid2, tagName, mpicom_join, source_id)
          if (ierr .ne. 0) then
             call shr_sys_abort(subname//' cannot receive element tag: '//trim(tagName))
          endif
       endif

       !---------------------------------------------------------------------------
       ! Free sender communication buffers to conserve memory
       ! Important for large-scale runs to avoid memory buildup
       !---------------------------------------------------------------------------
       if (mbAPPid1 .ge. 0) then
          ierr = iMOAB_FreeSenderBuffers(mbAPPid1, target_id)
          if (ierr .ne. 0) then
             call shr_sys_abort(subname//' cannot free sender buffers')
          endif
       endif

       !---------------------------------------------------------------------------
       ! Stop map exchange timer
       !---------------------------------------------------------------------------
       if (present(timer_map_exch)) then
          call t_drvstopf (trim(timer_map_exch))
       end if
    endif

    !---------------------------------------------------------------------------
    ! Exchange infodata (metadata) if provided
    ! Infodata contains runtime configuration and state information
    !---------------------------------------------------------------------------
    if (present(timer_infodata_exch)) then
       if (present(mpicom_barrier)) then
          call t_drvstartf (trim(timer_infodata_exch), barrier=mpicom_barrier)
       end if
    end if

    if (present(infodata) .and. present(infodata_string)) then
       if (direction == 'c2x') then  ! component to coupler
          if (comp%iamin_cplcompid) then
             call seq_infodata_exchange(infodata, comp%cplcompid, trim(infodata_string))
          end if
       else  ! x2c: coupler to component
          if (comp%iamin_cplallcompid) then
             call seq_infodata_exchange(infodata, comp%cplallcompid, trim(infodata_string))
          end if
       endif
    end if

    if (present(timer_infodata_exch)) then
       call t_drvstopf (trim(timer_infodata_exch))
    end if

    !---------------------------------------------------------------------------
    ! Stop component exchange timer
    !---------------------------------------------------------------------------
    if (present(timer_comp_exch)) then
       if (present(mpicom_barrier)) then
          call t_drvstopf (trim(timer_comp_exch), cplcom=.true.)
       end if
    end if

    !---------------------------------------------------------------------------
    ! Debug output: write mesh files when MOABDEBUG is defined
    ! Useful for visualizing data exchange and debugging coupling issues
    !---------------------------------------------------------------------------
#ifdef MOABDEBUG
    write(lnum,"(I0.2)") num_moab_exports
    if (seq_comm_iamroot(CPLID) ) then
       write(logunit,'(A)') subname//' '//comp%ntype//' at moab count '//trim(lnum)//' called in direction '//trim(direction)//' for fields '//trim(tagname)
    endif

    ! Write received mesh data to HDF5 file for visualization/debugging
    if (mbAPPid2 .ge. 0 ) then !  we are on receiving pes, for sure
      if (present(context_exch)) then
         outfile = comp%ntype//'_'//trim(context_exch)//'_'//trim(direction)//'_'//trim(lnum)//'.h5m'//C_NULL_CHAR
      else
         outfile = comp%ntype//'_'//trim(direction)//'_'//trim(lnum)//'.h5m'//C_NULL_CHAR
      endif
      wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(mbAPPid2, trim(outfile), trim(wopts))
      if (ierr .ne. 0) then
          call shr_sys_abort(subname//' cannot write file '// outfile)
       endif
    endif
#endif

  end subroutine component_exch_moab

end module cplcomp_exchange_mod
