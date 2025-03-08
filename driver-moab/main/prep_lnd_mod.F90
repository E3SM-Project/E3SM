module prep_lnd_mod

  use shr_kind_mod    , only: R8 => SHR_KIND_R8
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_kind_mod    , only: cxx => SHR_KIND_CXX
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_atm, num_inst_rof, num_inst_glc
  use seq_comm_mct    , only: num_inst_lnd, num_inst_frc
  use seq_comm_mct    , only: CPLID, LNDID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_comm_mct,     only: mlnid  ! iMOAB pid for land mesh on component pes
  use seq_comm_mct,     only: mhid     ! iMOAB id for atm instance
  use seq_comm_mct,     only: mphaid   ! iMOAB id for phys atm on atm pes
  use seq_comm_mct,     only: mhpgid   ! iMOAB id for atm pgx grid, on atm pes; created with se and gll grids
  use seq_comm_mct,     only: mblxid ! iMOAB id for land migrated mesh to coupler pes
  use seq_comm_mct,     only: mbrxid   !          iMOAB id of moab rof on coupler pes (FV now)
  use seq_comm_mct,     only: mbintxal ! iMOAB id for intx mesh between atm and lnd
  use seq_comm_mct,     only: mbintxrl ! iMOAB id for intx mesh between river and land
  use seq_comm_mct,     only: mb_rof_aream_computed  ! signal

  use seq_comm_mct,     only: mbaxid   ! iMOAB id for atm migrated mesh to coupler pes
  use seq_comm_mct,     only: atm_pg_active  ! whether the atm uses FV mesh or not ; made true if fv_nphys > 0
  ! use dimensions_mod,   only: np     ! for atmosphere
  use seq_comm_mct,     only: seq_comm_getinfo => seq_comm_setptrs
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: lnd, atm, rof, glc
  use map_glc2lnd_mod   , only: map_glc2lnd_ec
  use iso_c_binding
#ifdef  HAVE_MOAB
  use iMOAB , only: iMOAB_ComputeCommGraph, iMOAB_ComputeMeshIntersectionOnSphere, &
    iMOAB_ComputeScalarProjectionWeights, iMOAB_DefineTagStorage, iMOAB_RegisterApplication, &
    iMOAB_WriteMesh, iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorage
  use seq_comm_mct,     only : num_moab_exports
#endif

#ifdef MOABCOMP
  use component_type_mod, only:  compare_mct_av_moab_tag
#endif

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_lnd_init
  public :: prep_lnd_mrg
  ! moab version
  public :: prep_lnd_mrg_moab

  public :: prep_lnd_calc_a2x_lx
  public :: prep_lnd_calc_r2x_lx
  public :: prep_lnd_calc_g2x_lx
  public :: prep_lnd_calc_z2x_lx

  public :: prep_lnd_get_a2x_lx
  public :: prep_lnd_get_r2x_lx
  public :: prep_lnd_get_g2x_lx
  public :: prep_lnd_get_z2x_lx

  public :: prep_lnd_get_mapper_Sa2l
  public :: prep_lnd_get_mapper_Fa2l
  public :: prep_lnd_get_mapper_Fr2l
  public :: prep_lnd_get_mapper_Sg2l
  public :: prep_lnd_get_mapper_Fg2l

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_lnd_merge
  private :: prep_lnd_set_glc2lnd_fields

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sa2l           ! needed in ccsm_comp_mod.F90 (setting of aream)
  type(seq_map), pointer :: mapper_Fa2l           ! needed in ccsm_comp_mod.F90 (seq_domain_check)
  type(seq_map), pointer :: mapper_Fr2l           ! needed in seq_frac_mct.F90
  type(seq_map), pointer :: mapper_Sg2l           ! currently unused (all g2l mappings use the flux mapper)
  type(seq_map), pointer :: mapper_Fg2l

  ! attribute vectors
  type(mct_aVect), pointer :: a2x_lx(:) ! Atm export, lnd grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: r2x_lx(:) ! Rof export, lnd grid, lnd pes - allocated in lnd gc
  type(mct_aVect), pointer :: g2x_lx(:) ! Glc export, lnd grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: z2x_lx(:) ! Iac export, lnd grid, cpl pes - allocated in driver

  ! seq_comm_getData variables
  integer :: mpicom_CPLID                         ! MPI cpl communicator

  ! field names and lists, for fields that need to be treated specially
  character(len=*), parameter :: glc_frac_field = 'Sg_ice_covered'
  character(len=*), parameter :: glc_topo_field = 'Sg_topo'
  character(len=*), parameter :: glc_icemask_field = 'Sg_icemask'
  ! fields mapped from glc to lnd, NOT separated by elevation class
  character(CXX) :: glc2lnd_non_ec_fields
  ! other fields (besides frac_field and topo_field) that are mapped from glc to lnd,
  ! separated by elevation class
  character(CXX) :: glc2lnd_ec_extra_fields
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_lnd_init(infodata, atm_c2_lnd, rof_c2_lnd, glc_c2_lnd, iac_c2_lnd)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_lnd ! .true.  => atm to lnd coupling on
    logical                 , intent(in)    :: rof_c2_lnd ! .true.  => rof to lnd coupling on
    logical                 , intent(in)    :: glc_c2_lnd ! .true.  => glc to lnd coupling on
    logical                 , intent(in)    :: iac_c2_lnd ! .true.  => iac to lnd coupling on
    !
    ! Local Variables
    integer                  :: lsize_l
    integer                  :: eai, eri, egi
    logical                  :: samegrid_al   ! samegrid atm and land
    logical                  :: samegrid_lr   ! samegrid land and rof
    logical                  :: samegrid_lg   ! samegrid land and glc
    logical                  :: esmf_map_flag ! .true. => use esmf for mapping
    logical                  :: lnd_present   ! .true. => land is present
    logical                  :: iamroot_CPLID ! .true. => CPLID masterproc
    character(CL)            :: atm_gnam      ! atm grid
    character(CL)            :: lnd_gnam      ! lnd grid
    character(CL)            :: rof_gnam      ! rof grid
    character(CL)            :: glc_gnam      ! glc grid
    type(mct_avect), pointer :: l2x_lx
#ifdef HAVE_MOAB
   ! MOAB stuff
    integer                  :: ierr, idintx, rank
    character*32             :: appname, outfile, wopts, lnum
    character*32             :: dm1, dm2, dofnameS, dofnameT, wgtIdef
    integer                  :: orderS, orderT, volumetric, noConserve, validate, fInverseDistanceMap
    integer                  :: fNoBubble, monotonicity
    ! will do comm graph over coupler PES, in 2-hop strategy
    integer                  :: mpigrp_CPLID ! coupler pes group, used for comm graph  <-> atm-lnd, rof-lnd

    integer                  :: type1, type2 ! type for computing graph; should be the same type for ocean, 3 (FV)
    integer                  :: tagtype, numco, tagindex
    character(CXX)           :: tagName

    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
    integer  mlsize ! moab land size
    integer  nrflds  ! number of rof fields projected on land
    integer arrsize  ! for setting the r2x fields on land to 0
    integer ent_type ! for setting tags
    real (kind=R8) , allocatable :: tmparray (:) ! used to set the r2x fields to 0

#endif
    character(*), parameter  :: subname = '(prep_lnd_init)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         esmf_map_flag=esmf_map_flag,   &
         lnd_present=lnd_present,       &
         atm_gnam=atm_gnam,             &
         lnd_gnam=lnd_gnam,             &
         rof_gnam=rof_gnam,             &
         glc_gnam=glc_gnam)

    allocate(mapper_Sa2l)
    allocate(mapper_Fa2l)
    allocate(mapper_Fr2l)
    allocate(mapper_Sg2l)
    allocate(mapper_Fg2l)

    if (lnd_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       allocate(a2x_lx(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_lx(eai), rList=seq_flds_a2x_fields, lsize=lsize_l)
          call mct_aVect_zero(a2x_lx(eai))
       enddo
       allocate(r2x_lx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_aVect_init(r2x_lx(eri), rlist=seq_flds_r2x_fields, lsize=lsize_l)
          call mct_aVect_zero(r2x_lx(eri))
       end do
       allocate(g2x_lx(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(g2x_lx(egi), rList=seq_flds_x2l_fields_from_glc, lsize=lsize_l)
          call mct_aVect_zero(g2x_lx(egi))
       end do

       samegrid_al = .true.
       samegrid_lr = .true.
       samegrid_lg = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.
       if (trim(lnd_gnam) /= trim(glc_gnam)) samegrid_lg = .false.

       if (rof_c2_lnd) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fr2l'
          end if
          call seq_map_init_rcfile(mapper_Fr2l, rof(1), lnd(1), &
               'seq_maps.rc','rof2lnd_fmapname:','rof2lnd_fmaptype:',samegrid_lr, &
               string='mapper_Fr2l initialization',esmf_map=esmf_map_flag)
! symmetric of l2r, from prep_rof
#ifdef HAVE_MOAB
          ! Call moab intx only if land and river are init in moab
          if ((mbrxid .ge. 0) .and.  (mblxid .ge. 0)) then
            appname = "ROF_LND_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between rof and lnd mesh
            idintx = 100*rof(1)%cplcompid + lnd(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxrl)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering  rof lnd intx'
              call shr_sys_abort(subname//' ERROR in registering rof lnd intx')
            endif
            if (samegrid_lr)then
! the same mesh , lnd and rof use the same dofs, but restricted
               ! we do not compute intersection, so we will have to just send data from lnd to rof and viceversa, by GLOBAL_ID matching
               ! so we compute just a comm graph, between lnd and rof dofs, on the coupler; target is rof
               ! land is full mesh
               call seq_comm_getData(CPLID ,mpigrp=mpigrp_CPLID)
               type1 = 3; !  full mesh for lrofarnd now
               type2 = 3;  ! fv for target land
               ierr = iMOAB_ComputeCommGraph( mbrxid, mblxid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        rof(1)%cplcompid, lnd(1)%cplcompid)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing comm graph , rof-lnd'
                  call shr_sys_abort(subname//' ERROR in computing comm graph , rof-lnd')
               endif
               ! context for rearrange is target in this case
            if ( mapper_Fr2l%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fr2l%mbname) &
                             //' mapper_Fr2l'
                endif
            endif

               mapper_Fr2l%src_mbid = mbrxid
               mapper_Fr2l%tgt_mbid = mblxid
               mapper_Fr2l%src_context = rof(1)%cplcompid
               mapper_Fr2l%intx_context = lnd(1)%cplcompid
            else
              ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mbrxid, mblxid, mbintxrl)
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing   rof lnd intx'
                call shr_sys_abort(subname//' ERROR in computing  rof lnd intx')
              endif
              if (iamroot_CPLID) then
                write(logunit,*) 'iMOAB intersection between  rof and lnd with id:', idintx
              end if
              ! we also need to compute the comm graph for the second hop, from the rof on coupler to the
              ! rof for the intx rof-lnd context (coverage)
              !
              call seq_comm_getData(CPLID ,mpigrp=mpigrp_CPLID)
              type1 = 3 ! land is FV now on coupler side
              type2 = 3;

              ierr = iMOAB_ComputeCommGraph( mbrxid, mbintxrl, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                          rof(1)%cplcompid, idintx)
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing comm graph for second hop, lnd-rof'
                call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, lnd-rof')
              endif
              ! now take care of the mapper
            if ( mapper_Fr2l%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fr2l%mbname) &
                             //' mapper_Fr2l'
                endif
            endif
              mapper_Fr2l%src_mbid = mbrxid
              mapper_Fr2l%tgt_mbid = mblxid
              mapper_Fr2l%intx_mbid = mbintxrl
              mapper_Fr2l%src_context = rof(1)%cplcompid
              mapper_Fr2l%intx_context = idintx
              wgtIdef = 'scalar'//C_NULL_CHAR
              mapper_Fr2l%weight_identifier = wgtIdef
              mapper_Fr2l%mbname = 'mapper_Fr2l'

              ! because we will project fields from rof to lnd grid, we need to define
              !  the r2x fields to lnd grid on coupler side

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
              validate = 0 !! important
              fInverseDistanceMap = 0
              if (iamroot_CPLID) then
                write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxrl=', mbintxrl, ' wgtIdef=', wgtIdef, &
                    'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT)
              endif
              ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxrl, wgtIdef, &
                                                trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT) )

              ! signal that the aream for rof has been computed
              mb_rof_aream_computed = .true.
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing rl weights '
                call shr_sys_abort(subname//' ERROR in computing rl weights ')
              endif

#ifdef MOABDEBUG
              wopts = C_NULL_CHAR
              call shr_mpi_commrank( mpicom_CPLID, rank )
              if (rank .lt. 5) then
                write(lnum,"(I0.2)")rank !
                outfile = 'intx_rl_'//trim(lnum)// '.h5m' // C_NULL_CHAR
                ierr = iMOAB_WriteMesh(mbintxrl, outfile, wopts) ! write local intx file
                if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in writing intx rl file '
                  call shr_sys_abort(subname//' ERROR in writing intx rl file ')
                endif
              endif
#endif
            endif
            tagname = trim(seq_flds_r2x_fields)//C_NULL_CHAR
            tagtype = 1 ! dense
            numco = 1 !
            ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags for seq_flds_r2x_fields on lnd cpl'
               call shr_sys_abort(subname//' ERROR in  defining tags for seq_flds_r2x_fields on lnd cpl')
            endif

 ! find out the number of local elements in moab mesh land instance on coupler
            ierr  = iMOAB_GetMeshInfo ( mblxid, nvert, nvise, nbl, nsurf, nvisBC )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' cant get size of land mesh'
               call shr_sys_abort(subname//' ERROR in getting size of land mesh')
            endif
            ! land is now cell mesh on coupler side
            mlsize = nvise(1)
            ent_type = 1 ! cell
            ! set to 0 all fields that are projected from river
            nrflds = mct_aVect_nRattr(r2x_lx(1)) !  these are the numbers of fields in seq_flds_r2x_fields
            arrsize = nrflds*mlsize
            allocate (tmparray(arrsize)) ! mlsize is the size of local land
            ! do we need to zero out others or just river ?
            tmparray = 0._R8
            ierr = iMOAB_SetDoubleTagStorage(mblxid, tagname, arrsize , ent_type, tmparray)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' cant zero out r2x tags on land'
               call shr_sys_abort(subname//' cant zero out r2x tags on land')
            endif
            deallocate (tmparray)

         end if ! if ((mbrxid .ge. 0) .and.  (mblxid .ge. 0))
! endif HAVE_MOAB
#endif
       end if
       call shr_sys_flush(logunit)

       if (atm_c2_lnd) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2l'
          end if
          call seq_map_init_rcfile(mapper_Sa2l, atm(1), lnd(1), &
               'seq_maps.rc','atm2lnd_smapname:','atm2lnd_smaptype:',samegrid_al, &
               'mapper_Sa2l initialization',esmf_map_flag)
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fa2l'
          end if
          call seq_map_init_rcfile(mapper_Fa2l, atm(1), lnd(1), &
               'seq_maps.rc','atm2lnd_fmapname:','atm2lnd_fmaptype:',samegrid_al, &
               'mapper_Fa2l initialization',esmf_map_flag)
! similar to prep_atm_init, lnd and atm reversed
#ifdef HAVE_MOAB
          ! important change: do not compute intx at all between atm and land when we have samegrid_al
          ! we will use just a comm graph to send data from atm to land on coupler
          ! this is just a rearrange in a way
          if ((mbaxid .ge. 0) .and.  (mblxid .ge. 0) ) then
            appname = "ATM_LND_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between lnd and atm mesh
            idintx = 100*atm(1)%cplcompid + lnd(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxal)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering atm lnd intx '
              call shr_sys_abort(subname//' ERROR in registering atm lnd intx ')
            endif
            if ( mapper_Sa2l%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Sa2l%mbname) &
                             //' mapper_Sa2l'
                endif
            endif

            ! set up the scalar mapper context
            mapper_Sa2l%src_mbid = mbaxid
            mapper_Sa2l%tgt_mbid = mblxid
            mapper_Sa2l%intx_mbid = mbintxal
            mapper_Sa2l%src_context = atm(1)%cplcompid
            mapper_Sa2l%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Sa2l%weight_identifier = wgtIdef
            mapper_Sa2l%mbname = 'mapper_Sa2l'

            call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)
            if (.not. samegrid_al) then ! tri grid case
              if (iamroot_CPLID) then
                write(logunit,*) 'iMOAB intersection between atm and land with id:', idintx
              endif
              ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mbaxid, mblxid, mbintxal)
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing atm lnd intx'
                call shr_sys_abort(subname//' ERROR in computing atm lnd intx')
              endif
#ifdef MOABDEBUG
              ! write intx only if true intx file:
              wopts = C_NULL_CHAR
              call shr_mpi_commrank( mpicom_CPLID, rank )
                if (rank .lt. 3) then ! write only a few intx files
                write(lnum,"(I0.2)")rank !
                outfile = 'intx_al'//trim(lnum)// '.h5m' // C_NULL_CHAR
                ierr = iMOAB_WriteMesh(mbintxal, outfile, wopts) ! write local intx file
                if (ierr .ne. 0) then
                    write(logunit,*) subname,' error in writing intx file atm land  '
                    call shr_sys_abort(subname//' ERROR in writing intx file atm lnd')
                endif
              endif
#endif
              ! we also need to compute the comm graph for the second hop, from the atm on coupler to the
              ! lnd for the intx atm-lnd context (coverage)
              !
              if (atm_pg_active) then
                type1 = 3; !  fv for atm; cgll does not work anyway
              else
                type1 = 1 ! this projection works (cgll to fv), but reverse does not ( fv - cgll)
              endif
              type2 = 3; ! land is fv in this case (separate grid)

              ierr = iMOAB_ComputeCommGraph( mbaxid, mbintxal, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        atm(1)%cplcompid, idintx)
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing comm graph for second hop, atm-lnd'
                call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, atm-lnd')
              endif

              volumetric = 0 ! can be 1 only for FV->DGLL or FV->CGLL;

              if (atm_pg_active) then
                dm1 = "fv"//C_NULL_CHAR
                dofnameS="GLOBAL_ID"//C_NULL_CHAR
                orderS = 1 !  fv-fv
              else
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
                  write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxal=', mbintxal, ' wgtIdef=', wgtIdef, &
                      'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                                  fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                  noConserve, validate, &
                                                  trim(dofnameS), trim(dofnameT)
              endif
              ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxal, wgtIdef, &
                                                  trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                                  fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                  noConserve, validate, &
                                                  trim(dofnameS), trim(dofnameT) )
              if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing weights for atm-lnd   '
                  call shr_sys_abort(subname//' ERROR in computing weights for atm-lnd ')
              endif

            else  ! the same mesh , atm and lnd use the same dofs, but lnd is a subset of atm
                ! we do not compute intersection, so we will have to just send data from atm to land and viceversa, by GLOBAL_ID matching
                ! so we compute just a comm graph, between atm and lnd dofs, on the coupler; target is lnd
              ! land is point cloud in this case, type1 = 2

              if (atm_pg_active) then
                  type1 = 3; !  fv for atm; cgll does not work anyway
              else
                  type1 = 1 ! this projection works (cgll to fv), but reverse does not ( fv - cgll)
              endif
              type2 = 3;  ! FV mesh on coupler land
              ierr = iMOAB_ComputeCommGraph( mbaxid, mblxid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                      atm(1)%cplcompid, lnd(1)%cplcompid)
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing comm graph for second hop, atm-lnd'
                call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, atm-lnd')
              endif
              mapper_Sa2l%intx_context = lnd(1)%cplcompid

            endif ! if tri-grid

             ! use the same map for fluxes too
            if ( mapper_Fa2l%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fa2l%mbname) &
                             //' mapper_Fa2l'
                endif
            endif
            mapper_Fa2l%src_mbid = mbaxid
            mapper_Fa2l%tgt_mbid = mblxid
            mapper_Fa2l%intx_mbid = mbintxal
            mapper_Fa2l%src_context = atm(1)%cplcompid
            mapper_Fa2l%intx_context = mapper_Sa2l%intx_context
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fa2l%weight_identifier = wgtIdef
            mapper_Fa2l%mbname = 'mapper_Fa2l'


            ! in any case, we need to define the tags on landx from the phys atm seq_flds_a2x_fields
            tagtype = 1  ! dense, double
            numco = 1 !  one value per vertex / entity
            tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
            if ( ierr > 0) then
                call shr_sys_abort(subname//' fail to define seq_flds_a2x_fields for lnd x moab mesh ')
            endif

          endif    ! if ((mbaxid .ge. 0) .and.  (mblxid .ge. 0) ) then


#endif
       endif
       call shr_sys_flush(logunit)

       if (glc_c2_lnd) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sg2l'
          end if
          call seq_map_init_rcfile(mapper_Sg2l, glc(1), lnd(1), &
               'seq_maps.rc','glc2lnd_smapname:','glc2lnd_smaptype:',samegrid_lg, &
               'mapper_Sg2l initialization',esmf_map_flag)
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fg2l'
          end if
          call seq_map_init_rcfile(mapper_Fg2l, glc(1), lnd(1), &
               'seq_maps.rc','glc2lnd_fmapname:','glc2lnd_fmaptype:',samegrid_lg, &
               'mapper_Fg2l initialization',esmf_map_flag)

          call prep_lnd_set_glc2lnd_fields()
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_lnd_init

  !================================================================================================

  subroutine prep_lnd_set_glc2lnd_fields()

    !---------------------------------------------------------------
    ! Description
    ! Sets the module-level glc2lnd_non_ec_fields and glc2lnd_ec_extra_fields variables.
    !
    ! Local Variables
    character(len=CXX) :: temp_list

    character(*), parameter  :: subname = '(prep_lnd_set_glc2lnd_fields)'
    !---------------------------------------------------------------

    ! glc2lnd fields not separated by elevation class can be determined by finding fields
    ! that exist in both the g2x_to_lnd list and the x2l_from_glc list
    call shr_string_listIntersect(seq_flds_g2x_fields_to_lnd, &
         seq_flds_x2l_fields_from_glc, &
         glc2lnd_non_ec_fields)

    ! glc2lnd fields separated by elevation class are all fields not determined above.
    ! However, we also need to remove glc_frac_field and glc_topo_field from this list,
    ! because those are handled specially, so are not expected to be in this
    ! "extra_fields" list.
    !
    ! NOTE(wjs, 2015-04-24) I am going to the trouble of building this field list
    ! dynamically, rather than simply hard-coding the necessary fields (currently just
    ! 'Flgg_hflx'), so that new fields can be added in seq_flds_mod without needing to
    ! change any other code.
    call shr_string_listDiff(seq_flds_g2x_fields_to_lnd, &
         glc2lnd_non_ec_fields, &
         glc2lnd_ec_extra_fields)
    temp_list = glc2lnd_ec_extra_fields
    call shr_string_listDiff(temp_list, &
         glc_frac_field, &
         glc2lnd_ec_extra_fields)
    temp_list = glc2lnd_ec_extra_fields
    call shr_string_listDiff(temp_list, &
         glc_topo_field, &
         glc2lnd_ec_extra_fields)

  end subroutine prep_lnd_set_glc2lnd_fields

  !================================================================================================

  subroutine prep_lnd_mrg(infodata, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Prepare run phase, including running the merge
    !
    ! Arguments
    type(seq_infodata_type) , intent(in) :: infodata
    character(len=*)     , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eai, eri, egi, eli
    type(mct_aVect), pointer :: x2l_lx
    character(*), parameter  :: subname = '(prep_lnd_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       ! Use fortran mod to address ensembles in merge
       eai = mod((eli-1),num_inst_atm) + 1
       eri = mod((eli-1),num_inst_rof) + 1
       egi = mod((eli-1),num_inst_glc) + 1

       x2l_lx => component_get_x2c_cx(lnd(eli))  ! This is actually modifying x2l_lx
       call prep_lnd_merge( a2x_lx(eai), r2x_lx(eri), g2x_lx(egi), x2l_lx )
    enddo
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_lnd_mrg

! this does almost nothing now, except documenting
  subroutine prep_lnd_mrg_moab (infodata)
    type(seq_infodata_type) , intent(in) :: infodata


    type(mct_avect) , pointer   :: a2x_l  ! used just for indexing
    type(mct_avect) , pointer   :: r2x_l
    type(mct_avect) , pointer   :: g2x_l
    type(mct_avect) , pointer   :: x2l_l

    !--------------------------------------------------

    character(*), parameter  :: subname = '(prep_lnd_mrg_moab)'

    ! this routine does mostly nothing for moab, no fields are actually combined
    ! keep it here for documentation mostly
    ! Description
    ! Create input land state directly from atm, runoff and glc outputs
    !
    !-----------------------------------------------------------------------
    integer       :: nflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: r2x_sharedindices
    type(mct_aVect_sharedindices),save :: g2x_sharedindices
#ifdef MOABDEBUG
    integer :: ierr
    character*32             :: outfile, wopts, lnum
#endif
#ifdef MOABCOMP
    character(CXX)           :: tagname, mct_field
    real(R8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list, ent_type
    type(mct_string)    :: mctOStr  !
#endif
    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then
      a2x_l => a2x_lx(1)
      r2x_l => r2x_lx(1)
      g2x_l => g2x_lx(1)
      x2l_l => component_get_x2c_cx(lnd(1))
       nflds = mct_aVect_nRattr(x2l_l)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2l_l)
          mrgstr(i) = subname//'x2l%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(a2x_l, x2l_l, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(r2x_l, x2l_l, r2x_SharedIndices)
       call mct_aVect_setSharedIndices(g2x_l, x2l_l, g2x_SharedIndices)

       !--- document copy operations ---
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, a2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
       enddo
       do i=1,r2x_SharedIndices%shared_real%num_indices
          i1=r2x_SharedIndices%shared_real%aVindices1(i)
          o1=r2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, r2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = r2x%'//trim(field)
       enddo
       do i=1,g2x_SharedIndices%shared_real%num_indices
          i1=g2x_SharedIndices%shared_real%aVindices1(i)
          o1=g2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, g2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field)
       enddo
    endif

    ! call mct_aVect_copy(aVin=a2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    ! call mct_aVect_copy(aVin=r2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=r2x_SharedIndices)
    ! call mct_aVect_copy(aVin=g2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

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
  ! land does not do any merge for moab, all fields are directly projected, from atm, river, glacier
  ! compare_mct_av_moab_tag(comp, attrVect, field, imoabApp, tag_name, ent_type, difference)
    x2l_l => component_get_x2c_cx(lnd(1))
    ! loop over all fields in seq_flds_x2l_fields
    call mct_list_init(temp_list ,seq_flds_x2l_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! cell for land now, it is a full mesh
    if (iamroot) print *, num_moab_exports, trim(seq_flds_x2l_fields)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call compare_mct_av_moab_tag(lnd(1), x2l_l, mct_field,  mblxid, tagname, ent_type, difference, first_time)
    enddo
    call mct_list_clean(temp_list)
#endif

    first_time = .false.

#ifdef MOABDEBUG
    if (mblxid .ge. 0 ) then !  we are on coupler pes, for sure
       write(lnum,"(I0.2)")num_moab_exports
       outfile = 'LndCplAftMm'//trim(lnum)//'.h5m'//C_NULL_CHAR
       wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
       ierr = iMOAB_WriteMesh(mblxid, trim(outfile), trim(wopts))
    endif

#endif

  end subroutine prep_lnd_mrg_moab
  !================================================================================================

  subroutine prep_lnd_merge( a2x_l, r2x_l, g2x_l, x2l_l )
    !---------------------------------------------------------------
    ! Description
    ! Create input land state directly from atm, runoff and glc outputs
    !
    ! Arguments
    type(mct_aVect), intent(in)     :: a2x_l
    type(mct_aVect), intent(in)     :: r2x_l
    type(mct_aVect), intent(in)     :: g2x_l
    type(mct_aVect), intent(inout)  :: x2l_l
    !-----------------------------------------------------------------------
    integer       :: nflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: r2x_sharedindices
    type(mct_aVect_sharedindices),save :: g2x_sharedindices
    character(*), parameter   :: subname = '(prep_lnd_merge) '

    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then
       nflds = mct_aVect_nRattr(x2l_l)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2l_l)
          mrgstr(i) = subname//'x2l%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(a2x_l, x2l_l, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(r2x_l, x2l_l, r2x_SharedIndices)
       call mct_aVect_setSharedIndices(g2x_l, x2l_l, g2x_SharedIndices)

       !--- document copy operations ---
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, a2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
       enddo
       do i=1,r2x_SharedIndices%shared_real%num_indices
          i1=r2x_SharedIndices%shared_real%aVindices1(i)
          o1=r2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, r2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = r2x%'//trim(field)
       enddo
       do i=1,g2x_SharedIndices%shared_real%num_indices
          i1=g2x_SharedIndices%shared_real%aVindices1(i)
          o1=g2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, g2x_l)
          mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field)
       enddo
    endif

    call mct_aVect_copy(aVin=a2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    call mct_aVect_copy(aVin=r2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=r2x_SharedIndices)
    call mct_aVect_copy(aVin=g2x_l, aVout=x2l_l, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

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

  end subroutine prep_lnd_merge

  !================================================================================================

  subroutine prep_lnd_calc_a2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create  a2x_lx (note that a2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    type(mct_aVect), pointer :: a2x_ax
    character(*), parameter  :: subname = '(prep_lnd_calc_a2x_lx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       a2x_ax => component_get_c2x_cx(atm(eai))
       call seq_map_map(mapper_Fa2l, a2x_ax, a2x_lx(eai), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_a2x_lx

  !================================================================================================

  subroutine prep_lnd_calc_r2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create r2x_lx (note that r2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri
    type(mct_aVect) , pointer :: r2x_rx
    character(*), parameter :: subname = '(prep_lnd_calc_r2x_lx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       r2x_rx => component_get_c2x_cx(rof(eri))

       ! Note that one of these fields (a volr field) is remapped from rof -> lnd in
       ! map_lnd2rof_irrig_mod, because it is needed as a normalization term. So, if the
       ! details of this mapping call are changed in the future, it's possible that the
       ! equivalent r2l mapping in map_lnd2rof_irrig_mod should be changed to keep the two
       ! equivalent.
       call seq_map_map(mapper_Fr2l, r2x_rx, r2x_lx(eri), &
            fldlist=seq_flds_r2x_fluxes, norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_r2x_lx

  !================================================================================================

  subroutine prep_lnd_calc_g2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_lx (note that g2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_aVect), pointer :: g2x_gx
    character(*), parameter :: subname = '(prep_lnd_calc_g2x_lx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       g2x_gx => component_get_c2x_cx(glc(egi))

       ! Map fields that are NOT separated by elevation class on the land grid
       !
       ! These are mapped using a simple area-conservative remapping. (Note that we use
       ! the flux mapper even though these contain states, because we need these icemask
       ! fields to be mapped conservatively.)
       !
       ! Note that this mapping is redone for Sg_icemask in prep_glc_mod:
       ! prep_glc_map_qice_conservative_lnd2glc. If we ever change this mapping (e.g.,
       ! changing norm to .false.), then we should change the mapping there, too.
       !
       ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but
       ! this requires some more thought
       call seq_map_map(mapper_Fg2l, g2x_gx, g2x_lx(egi), &
            fldlist = glc2lnd_non_ec_fields, norm=.true.)

       ! Map fields that are separated by elevation class on the land grid
       call map_glc2lnd_ec( &
            g2x_g = g2x_gx, &
            frac_field = glc_frac_field, &
            topo_field = glc_topo_field, &
            icemask_field = glc_icemask_field, &
            extra_fields = glc2lnd_ec_extra_fields, &
            mapper = mapper_Fg2l, &
            g2x_l = g2x_lx(egi))
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_lnd_calc_g2x_lx

  !================================================================================================

  subroutine prep_lnd_calc_z2x_lx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create z2x_lx (note that z2x_lx is a local module variable)
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_aVect), pointer :: z2x_gx
    character(*), parameter :: subname = '(prep_lnd_calc_z2x_lx)'
    !---------------------------------------------------------------

    ! Stub

  end subroutine prep_lnd_calc_z2x_lx

  !================================================================================================

  function prep_lnd_get_a2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_a2x_lx(:)
    prep_lnd_get_a2x_lx => a2x_lx(:)
  end function prep_lnd_get_a2x_lx

  function prep_lnd_get_r2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_r2x_lx(:)
    prep_lnd_get_r2x_lx => r2x_lx(:)
  end function prep_lnd_get_r2x_lx

  function prep_lnd_get_g2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_g2x_lx(:)
    prep_lnd_get_g2x_lx => g2x_lx(:)
  end function prep_lnd_get_g2x_lx

  function prep_lnd_get_z2x_lx()
    type(mct_aVect), pointer :: prep_lnd_get_z2x_lx(:)
    prep_lnd_get_z2x_lx => z2x_lx(:)
  end function prep_lnd_get_z2x_lx

  function prep_lnd_get_mapper_Sa2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Sa2l
    prep_lnd_get_mapper_Sa2l => mapper_Sa2l
  end function prep_lnd_get_mapper_Sa2l

  function prep_lnd_get_mapper_Fa2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Fa2l
    prep_lnd_get_mapper_Fa2l => mapper_Fa2l
  end function prep_lnd_get_mapper_Fa2l

  function prep_lnd_get_mapper_Fr2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Fr2l
    prep_lnd_get_mapper_Fr2l => mapper_Fr2l
  end function prep_lnd_get_mapper_Fr2l

  function prep_lnd_get_mapper_Sg2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Sg2l
    prep_lnd_get_mapper_Sg2l => mapper_Sg2l
  end function prep_lnd_get_mapper_Sg2l

  function prep_lnd_get_mapper_Fg2l()
    type(seq_map), pointer :: prep_lnd_get_mapper_Fg2l
    prep_lnd_get_mapper_Fg2l => mapper_Fg2l
  end function prep_lnd_get_mapper_Fg2l

end module prep_lnd_mod
