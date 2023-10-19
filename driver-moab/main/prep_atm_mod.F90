module prep_atm_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: CX => shr_kind_CX, CXX => shr_kind_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_atm, num_inst_ocn, num_inst_ice, num_inst_lnd, num_inst_xao, &
       num_inst_frc, num_inst_max, CPLID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: atm, lnd, ocn, ice

  use shr_mpi_mod,  only:  shr_mpi_commrank
  use seq_comm_mct, only : mbaxid   ! iMOAB id for atm migrated mesh to coupler pes
  use seq_comm_mct, only : mphaid   ! iMOAB id for phys atm on atm pes
  use seq_comm_mct, only : mboxid   ! iMOAB id for mpas ocean migrated mesh to coupler pes
  use seq_comm_mct, only : mbofxid   ! iMOAB id for mpas ocean migrated mesh to coupler pes, just for xao flux calculations
  use seq_comm_mct, only : mbintxoa ! iMOAB id for intx mesh between ocean and atmosphere
  use seq_comm_mct, only : mbintxao ! iMOAB id for intx mesh between atm and ocean
  use seq_comm_mct, only : mbixid   ! iMOAB id for mpas ice migrated mesh to coupler pes
  use seq_comm_mct, only : mbintxia  ! iMOAB id for intx mesh between ice and atm
  use seq_comm_mct, only : mhid     ! iMOAB id for atm instance
  use seq_comm_mct, only : mhpgid   ! iMOAB id for atm pgx grid, on atm pes; created with se and gll grids
  use seq_comm_mct, only : atm_pg_active  ! whether the atm uses FV mesh or not ; made true if fv_nphys > 0
  use seq_comm_mct, only : mblxid   ! iMOAB id for land migrated to coupler pes !! old name : mlnxid
  use seq_comm_mct, only : mbintxla ! iMOAB id for intx mesh between land and atmosphere
  use seq_comm_mct, only : seq_comm_getinfo => seq_comm_setptrs
  use seq_comm_mct, only : num_moab_exports

  use dimensions_mod, only : np     ! for atmosphere
#ifdef MOABCOMP
  use component_type_mod, only:  compare_mct_av_moab_tag
#endif

  use iso_c_binding


  implicit none
  save
  PRIVATE

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_atm_init
  public :: prep_atm_mrg
  public :: prep_atm_mrg_moab

  public :: prep_atm_get_l2x_ax
  public :: prep_atm_get_i2x_ax
  public :: prep_atm_get_o2x_ax
  public :: prep_atm_get_z2x_ax

  public :: prep_atm_calc_l2x_ax
  public :: prep_atm_calc_i2x_ax
  public :: prep_atm_calc_o2x_ax
  public :: prep_atm_calc_z2x_ax

  public :: prep_atm_get_mapper_So2a
  public :: prep_atm_get_mapper_Fo2a
  public :: prep_atm_get_mapper_Sof2a
  public :: prep_atm_get_mapper_Fof2a
  public :: prep_atm_get_mapper_Sl2a
  public :: prep_atm_get_mapper_Fl2a
  public :: prep_atm_get_mapper_Si2a
  public :: prep_atm_get_mapper_Fi2a


  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_atm_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_So2a
  type(seq_map), pointer :: mapper_Sof2a          ! for moab fluxes
  type(seq_map), pointer :: mapper_Sl2a
  type(seq_map), pointer :: mapper_Si2a
  type(seq_map), pointer :: mapper_Fo2a           ! needed for seq_frac_init
  type(seq_map), pointer :: mapper_Fof2a
  type(seq_map), pointer :: mapper_Fl2a           ! needed for seq_frac_init
  type(seq_map), pointer :: mapper_Fi2a           ! needed for seq_frac_init

  ! attribute vectors
  type(mct_aVect), pointer :: l2x_ax(:)   ! Lnd export, atm grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: i2x_ax(:)   ! Ice export, atm grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: o2x_ax(:)   ! Ocn export, atm grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: z2x_ax(:)   ! Iac export, atm grid, cpl pes - allocated in driver

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
  logical :: iamroot_CPLID ! .true. => CPLID masterproc

#ifdef HAVE_MOAB
  real (kind=r8) , allocatable, private :: fractions_am (:,:) ! will retrieve the fractions from atm, and use them
  !  they were init with
  ! character(*),parameter :: fraclist_a = 'afrac:ifrac:ofrac:ifrad:ofrad' in moab, on the fractions
  real (kind=r8) , allocatable, private :: x2a_am (:,:)
  real (kind=r8) , allocatable, private :: l2x_am (:,:)
  real (kind=r8) , allocatable, private :: i2x_am (:,:)
  real (kind=r8) , allocatable, private :: o2x_am (:,:)
  !real (kind=r8) , allocatable, private :: z2x_am (:,:)
  real (kind=r8) , allocatable, private :: xao_am (:,:)  ! ?
#endif
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_atm_init(infodata, ocn_c2_atm, ice_c2_atm, lnd_c2_atm, iac_c2_atm)

   use iMOAB, only: iMOAB_ComputeMeshIntersectionOnSphere, iMOAB_RegisterApplication, &
   iMOAB_WriteMesh , iMOAB_ComputeCommGraph, iMOAB_ComputeScalarProjectionWeights, &
   iMOAB_DefineTagStorage
   !---------------------------------------------------------------
   ! Description
   ! Initialize module attribute vectors and  mappers
   !
   ! Arguments
   type (seq_infodata_type) , intent(inout) :: infodata
   logical                  , intent(in)    :: ocn_c2_atm ! .true.  => ocn to atm coupling on
   logical                  , intent(in)    :: ice_c2_atm ! .true.  => ice to atm coupling on
   logical                  , intent(in)    :: lnd_c2_atm ! .true.  => lnd to atm coupling on
   logical                  , intent(in)    :: iac_c2_atm ! .true.  => iac to atm coupling on
   !
   ! Local Variables
   integer                          :: lsize_a
   integer                          :: eli, eii, emi
   logical                          :: samegrid_ao    ! samegrid atm and ocean
   logical                          :: samegrid_al    ! samegrid atm and land
   logical                          :: esmf_map_flag  ! .true. => use esmf for mapping
   logical                          :: atm_present    ! .true.  => atm is present
   logical                          :: ocn_present    ! .true.  => ocn is present
   logical                          :: ice_present    ! .true.  => ice is present
   logical                          :: lnd_present    ! .true.  => lnd is prsent
   character(CL)                    :: ocn_gnam       ! ocn grid
   character(CL)                    :: atm_gnam       ! atm grid
   character(CL)                    :: lnd_gnam       ! lnd grid
   type(mct_avect), pointer         :: a2x_ax
   character(*), parameter          :: subname = '(prep_atm_init)'
   character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"

   ! MOAB stuff
   integer                  :: ierr, idintx, rank
   character*32             :: appname, outfile, wopts, lnum
   character*32             :: dm1, dm2, dofnameS, dofnameT, wgtIdef
   integer                  :: orderS, orderT, volumetric, noConserve, validate, fInverseDistanceMap
   integer                  :: fNoBubble, monotonicity
! will do comm graph over coupler PES, in 2-hop strategy
   integer                  :: mpigrp_CPLID ! coupler pes group, used for comm graph phys <-> atm-ocn

   integer                  :: type1, type2 ! type for computing graph; should be the same type for ocean, 3 (FV)
   integer                  :: tagtype, numco, tagindex
   character(CXX)           :: tagName
   integer                  :: context_id ! we will use a special context for the extra flux ocean instance
   logical                  :: no_match ! used to force a new mapper

   !---------------------------------------------------------------


   call seq_infodata_getData(infodata, &
      atm_present=atm_present,       &
      ocn_present=ocn_present,       &
      ice_present=ice_present,       &
      lnd_present=lnd_present,       &
      atm_gnam=atm_gnam,             &
      ocn_gnam=ocn_gnam,             &
      lnd_gnam=lnd_gnam,             &
      esmf_map_flag=esmf_map_flag)

   allocate(mapper_So2a)
   allocate(mapper_Sof2a)
   allocate(mapper_Sl2a)
   allocate(mapper_Si2a)
   allocate(mapper_Fo2a)
   allocate(mapper_Fof2a)
   allocate(mapper_Fl2a)
   allocate(mapper_Fi2a)

   if (atm_present) then

      call seq_comm_getData(CPLID, &
         mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

      a2x_ax => component_get_c2x_cx(atm(1))
      lsize_a = mct_aVect_lsize(a2x_ax)

      allocate(l2x_ax(num_inst_lnd))
      do eli = 1,num_inst_lnd
         call mct_aVect_init(l2x_ax(eli), rList=seq_flds_l2x_fields, lsize=lsize_a)
         call mct_aVect_zero(l2x_ax(eli))
      end do
      allocate(o2x_ax(num_inst_max))
      do emi = 1,num_inst_max
         call mct_aVect_init(o2x_ax(emi), rList=seq_flds_o2x_fields, lsize=lsize_a)
         call mct_aVect_zero(o2x_ax(emi))
      enddo
      allocate(i2x_ax(num_inst_ice))
      do eii = 1,num_inst_ice
         call mct_aVect_init(i2x_ax(eii), rList=seq_flds_i2x_fields, lsize=lsize_a)
         call mct_aVect_zero(i2x_ax(eii))
      enddo

      samegrid_al = .true.
      samegrid_ao = .true.
      if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
      if (trim(atm_gnam) /= trim(ocn_gnam)) samegrid_ao = .false.

      if (ocn_c2_atm) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_So2a'
         endif
         call seq_map_init_rcfile(mapper_So2a, ocn(1), atm(1), &
            'seq_maps.rc','ocn2atm_smapname:','ocn2atm_smaptype:',samegrid_ao, &
            'mapper_So2a initialization',esmf_map_flag)

         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Sof2a'
         endif
         no_match = .true. ! force to create a new mapper object
         call seq_map_init_rcfile(mapper_Sof2a, ocn(1), atm(1), &
            'seq_maps.rc','ocn2atm_smapname:','ocn2atm_smaptype:',samegrid_ao, &
            'mapper_Sof2a initialization',esmf_map_flag, no_match)

#ifdef HAVE_MOAB
         ! Call moab intx only if atm and ocn are init in moab
         if ((mbaxid .ge. 0) .and.  (mboxid .ge. 0)) then
            appname = "OCN_ATM_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between atm and ocn mesh
            idintx = 100*ocn(1)%cplcompid + atm(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxoa)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering atm ocn intx'
              call shr_sys_abort(subname//' ERROR in registering atm ocn intx')
            endif
            ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mboxid, mbaxid, mbintxoa)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in computing ocn atm intx'
              call shr_sys_abort(subname//' ERROR in computing ocn atm intx')
            endif
            if (iamroot_CPLID) then
              write(logunit,*) 'iMOAB intersection between ocean and atm with id:', idintx
            endif


            ! we also need to compute the comm graph for the second hop, from the ocn on coupler to the
            ! ocean for the intx ocean-atm context (coverage)
            !
            call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)
            type1 = 3; !  fv for ocean and atm; fv-cgll does not work anyway
            type2 = 3;
            ! ierr      = iMOAB_ComputeCommGraph( mboxid, mbintxoa, &mpicom_CPLID, &mpigrp_CPLID, &mpigrp_CPLID, &type1, &type2,
            !                              &ocn_id, &idintx)
            ierr = iMOAB_ComputeCommGraph( mboxid, mbintxoa, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        ocn(1)%cplcompid, idintx)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing comm graph for second hop, ocn-atm'
               call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, ocn-atm')
            endif

            ! now take care of the mapper
            if ( mapper_So2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_So2a%mbname) &
                             //' mapper_So2a'
                endif
            endif
            mapper_So2a%src_mbid = mboxid
            mapper_So2a%tgt_mbid = mbaxid !
            mapper_So2a%intx_mbid = mbintxoa
            mapper_So2a%src_context = ocn(1)%cplcompid
            mapper_So2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_So2a%weight_identifier = wgtIdef
            mapper_So2a%mbname = 'mapper_So2a'

            ! because we will project fields from ocean to atm phys grid, we need to define
            ! ocean o2x fields to atm phys grid (or atm spectral ext ) on coupler side

            if (atm_pg_active) then
               tagname = trim(seq_flds_o2x_fields)//C_NULL_CHAR
               tagtype = 1 ! dense
               numco = 1 !
               ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in defining tags for seq_flds_o2x_fields'
                  call shr_sys_abort(subname//' ERROR in coin defining tags for seq_flds_o2x_fields')
               endif
            else ! spectral case, fix later TODO
               numco = np*np !
            endif !

            volumetric = 0 ! can be 1 only for FV->DGLL or FV->CGLL;

            if (atm_pg_active) then
              dm2 = "fv"//C_NULL_CHAR
              dofnameT="GLOBAL_ID"//C_NULL_CHAR
              orderT = 1 !  fv-fv
            else
              dm2 = "cgll"//C_NULL_CHAR
              dofnameT="GLOBAL_DOFS"//C_NULL_CHAR
              orderT = np !  it should be 4
            endif
            dm1 = "fv"//C_NULL_CHAR
            dofnameS="GLOBAL_ID"//C_NULL_CHAR
            orderS = 1  !  not much arguing
            fNoBubble = 1
            monotonicity = 0 !
            noConserve = 0
            validate = 0 ! less verbose
            fInverseDistanceMap = 0
            if (iamroot_CPLID) then
               write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxoa=', mbintxoa, ' wgtIdef=', wgtIdef, &
                   'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameS), trim(dofnameT)
            endif
            ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxoa, wgtIdef, &
                                               trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameS), trim(dofnameT) )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in iMOAB_ComputeScalarProjectionWeights ocn atm   '
               call shr_sys_abort(subname//' ERROR in iMOAB_ComputeScalarProjectionWeights ocn atm ')
            endif


#ifdef MOABDEBUG
            wopts = C_NULL_CHAR
            call shr_mpi_commrank( mpicom_CPLID, rank )
            if (rank .lt. 3) then
              write(lnum,"(I0.2)")rank !
              outfile = 'intx_oa_'//trim(lnum)// '.h5m' // C_NULL_CHAR
              ierr = iMOAB_WriteMesh(mbintxoa, outfile, wopts) ! write local intx file
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in writing intx file '
                call shr_sys_abort(subname//' ERROR in writing intx file ')
              endif
            endif
! endif for MOABDEBUG
#endif
         endif ! if ((mbaxid .ge. 0) .and.  (mboxid .ge. 0)) then

! FLUX make the app and mapper for the a2o flux mappings
         if ((mbaxid .ge. 0) .and.  (mbofxid .ge. 0)) then
            ! we also need to compute the comm graph for the second hop, from the ocn on coupler to the
            ! ocean for the intx ocean-atm context (coverage)
            !
            call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)
            type1 = 3; !  fv for ocean and atm; fv-cgll does not work anyway
            type2 = 3;
            ! we ideintified the app mbofxid with !id_join = id_join + 1000! kind of random
            ! line 1267 in cplcomp_exchange_mod.F90
            context_id = ocn(1)%cplcompid + 1000
            ierr = iMOAB_ComputeCommGraph( mbofxid, mbintxoa, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        context_id, idintx)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing comm graph for second hop, ocnf -atm'
               call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, ocnf-atm')
            endif

            if ( mapper_Sof2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Sof2a%mbname) &
                             //' mapper_Sof2a'
                endif
            endif
            mapper_Sof2a%src_mbid = mbofxid
            mapper_Sof2a%tgt_mbid = mbaxid
            mapper_Sof2a%intx_mbid = mbintxoa
            mapper_Sof2a%src_context = context_id
            mapper_Sof2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Sof2a%weight_identifier = wgtIdef
            mapper_Sof2a%mbname = 'mapper_Sof2a'
         endif

! endif for HAVE_MOAB
#endif

      endif  ! if (ocn_c2_atm) then

      ! needed for domain checking
      if (ocn_present) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Fo2a'
         endif
         call seq_map_init_rcfile(mapper_Fo2a, ocn(1), atm(1), &
            'seq_maps.rc','ocn2atm_fmapname:','ocn2atm_fmaptype:',samegrid_ao, &
            'mapper_Fo2a initialization',esmf_map_flag)

         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Fof2a'
         endif
         no_match = .true. ! force to create a new mapper object
         call seq_map_init_rcfile(mapper_Fof2a, ocn(1), atm(1), &
            'seq_maps.rc','ocn2atm_fmapname:','ocn2atm_fmaptype:',samegrid_ao, &
            'mapper_Fof2a initialization',esmf_map_flag, no_match)

! copy mapper_So2a , maybe change the matrix ? still based on intersection ?
#ifdef HAVE_MOAB
         if ((mbaxid .ge. 0) .and.  (mboxid .ge. 0)) then
            ! now take care of the mapper
            if ( mapper_Fo2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fo2a%mbname) &
                             //' mapper_Fo2a'
                endif
            endif
            mapper_Fo2a%src_mbid = mboxid
            mapper_Fo2a%tgt_mbid = mbaxid
            mapper_Fo2a%intx_mbid = mbintxoa
            mapper_Fo2a%src_context = ocn(1)%cplcompid
            mapper_Fo2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fo2a%weight_identifier = wgtIdef
            mapper_Fo2a%mbname = 'mapper_Fo2a'
         endif
         if ((mbaxid .ge. 0) .and.  (mbofxid .ge. 0)) then
            if ( mapper_Fof2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fof2a%mbname) &
                             //' mapper_Fof2a'
                endif
            endif
            mapper_Fof2a%src_mbid = mbofxid
            mapper_Fof2a%tgt_mbid = mbaxid
            mapper_Fof2a%intx_mbid = mbintxoa
            mapper_Fof2a%src_context = ocn(1)%cplcompid
            mapper_Fof2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fof2a%weight_identifier = wgtIdef
            mapper_Fof2a%mbname = 'mapper_Fof2a'
         endif
! endif for HAVE_MOAB
#endif

      endif ! endif (ocn_present) then
      call shr_sys_flush(logunit)

      if (ice_c2_atm) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Si2a'
         endif
         no_match = .true. ! force to create a new mapper object
         ! otherwise it may find ocean map, and this will not work on ice vars
         call seq_map_init_rcfile(mapper_Si2a, ice(1), atm(1), &
            'seq_maps.rc','ice2atm_smapname:','ice2atm_smaptype:',samegrid_ao, &
            'mapper_Si2a initialization',esmf_map_flag, no_match)
            ! similar to ocn-atm mapping, do ice 2 atm mapping / set up

#ifdef HAVE_MOAB
         ! Call moab intx only if atm and ice are init in moab coupler
         if ((mbaxid .ge. 0) .and.  (mbixid .ge. 0)) then
            appname = "ICE_ATM_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between ice and atm mesh
            idintx = 100*ice(1)%cplcompid + atm(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxia)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering ice atm intx'
              call shr_sys_abort(subname//' ERROR in registering ice atm intx')
            endif
            ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mbixid, mbaxid, mbintxia)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in computing ice atm intx'
              call shr_sys_abort(subname//' ERROR in computing ice atm intx')
            endif
            if (iamroot_CPLID) then
              write(logunit,*) 'iMOAB intersection between ice and atm with id:', idintx
            endif


            ! we also need to compute the comm graph for the second hop, from the ice on coupler to the
            ! ice for the intx ice-atm context (coverage)
            !
            call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)
            type1 = 3; !  fv for ice and atm; fv-cgll does not work anyway
            type2 = 3;
            ! ierr      = iMOAB_ComputeCommGraph( mboxid, mbintxoa, &mpicom_CPLID, &mpigrp_CPLID, &mpigrp_CPLID, &type1, &type2,
            !                              &ocn_id, &idintx)
            ierr = iMOAB_ComputeCommGraph( mbixid, mbintxia, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        ice(1)%cplcompid, idintx)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing comm graph for second hop, ice-atm'
               call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, ice-atm')
            endif
            ! now take care of the mapper
            if ( mapper_Si2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Si2a%mbname) &
                             //' mapper_Si2a'
                endif
            endif
            mapper_Si2a%src_mbid = mbixid
            mapper_Si2a%tgt_mbid = mbaxid
            mapper_Si2a%intx_mbid = mbintxia
            mapper_Si2a%src_context = ice(1)%cplcompid
            mapper_Si2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Si2a%weight_identifier = wgtIdef
            mapper_Si2a%mbname = 'mapper_Si2a'
            ! because we will project fields from ocean to atm phys grid, we need to define
            ! ice i2x fields to atm phys grid (or atm spectral ext ) on coupler side
            if (atm_pg_active) then
               tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
               tagtype = 1 ! dense
               numco = 1 !
               ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in defining tags for seq_flds_i2x_fields'
                  call shr_sys_abort(subname//' ERROR in coin defining tags for seq_flds_i2x_fields')
               endif
            else ! spectral case, TODO
               tagtype = 1 ! dense
            endif

            volumetric = 0 ! can be 1 only for FV->DGLL or FV->CGLL;

            if (atm_pg_active) then
              dm2 = "fv"//C_NULL_CHAR
              dofnameT="GLOBAL_ID"//C_NULL_CHAR
              orderT = 1 !  fv-fv
            else
              dm2 = "cgll"//C_NULL_CHAR
              dofnameT="GLOBAL_DOFS"//C_NULL_CHAR
              orderT = np !  it should be 4
            endif
            dm1 = "fv"//C_NULL_CHAR
            dofnameS="GLOBAL_ID"//C_NULL_CHAR
            orderS = 1  !  not much arguing
            fNoBubble = 1
            monotonicity = 0 !
            noConserve = 0
            validate = 0 ! less verbose
            fInverseDistanceMap = 0
            if (iamroot_CPLID) then
               write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxia=', mbintxia, ' wgtIdef=', wgtIdef, &
                   'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameS), trim(dofnameT)
            endif
            ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxia, wgtIdef, &
                                               trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                               fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                               noConserve, validate, &
                                               trim(dofnameS), trim(dofnameT) )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in iMOAB_ComputeScalarProjectionWeights ice atm '
               call shr_sys_abort(subname//' error in iMOAB_ComputeScalarProjectionWeights ice atm ')
            endif


#ifdef MOABDEBUG
            wopts = C_NULL_CHAR
            call shr_mpi_commrank( mpicom_CPLID, rank )
            if (rank .lt. 3) then
               write(lnum,"(I0.2)")rank !
               outfile = 'intx_ia_'//trim(lnum)// '.h5m' // C_NULL_CHAR
               ierr = iMOAB_WriteMesh(mbintxia, outfile, wopts) ! write local intx file
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in writing intx file ice-atm '
                  call shr_sys_abort(subname//' ERROR in writing intx file ice-atm ')
               endif
            endif
! endif for MOABDEBUG
#endif
         endif ! if ((mbaxid .ge. 0) .and.  (mbixid .ge. 0)) then
! endif for HAVE_MOAB
#endif

      endif ! if (ice_c2_atm) then

      ! needed for domain checking
      if (ice_present) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Fi2a'
         endif
         no_match = .true. ! force a different map, we do not want to match to ocean
         call seq_map_init_rcfile(mapper_Fi2a, ice(1), atm(1), &
            'seq_maps.rc','ice2atm_fmapname:','ice2atm_fmaptype:',samegrid_ao, &
            'mapper_Fi2a initialization',esmf_map_flag, no_match)

#ifdef HAVE_MOAB
           ! now take care of the mapper for MOAB
            if ( mapper_Fi2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fi2a%mbname) &
                             //' mapper_Fi2a'
                endif
            endif
                
            mapper_Fi2a%src_mbid = mbixid
            mapper_Fi2a%tgt_mbid = mbaxid
            mapper_Fi2a%intx_mbid = mbintxia
            mapper_Fi2a%src_context = ice(1)%cplcompid
            mapper_Fi2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fi2a%weight_identifier = wgtIdef
            mapper_Fi2a%mbname = 'mapper_Fi2a'
#endif
      endif !  if (ice_present) then
      call shr_sys_flush(logunit)

      ! needed for domain checking
      if (lnd_present) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Fl2a'
         endif
         call seq_map_init_rcfile(mapper_Fl2a, lnd(1), atm(1), &
            'seq_maps.rc','lnd2atm_fmapname:','lnd2atm_fmaptype:',samegrid_al, &
            'mapper_Fl2a initialization',esmf_map_flag)

#ifdef HAVE_MOAB
         ! important change: do not compute intx at all between atm and land when we have samegrid_al
         ! we will use just a comm graph to send data from phys grid to land on coupler
         ! this is just a rearrange in a way
         if ((mbaxid .ge. 0) .and.  (mblxid .ge. 0) ) then
            appname = "LND_ATM_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between lnd and atm mesh
            idintx = 100*lnd(1)%cplcompid + atm(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxla)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering lnd atm intx '
              call shr_sys_abort(subname//' ERROR in registering lnd atm intx ')
            endif
            if ( mapper_Fl2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fl2a%mbname) &
                             //' mapper_Fl2a'
                endif
            endif
            mapper_Fl2a%src_mbid = mblxid
            mapper_Fl2a%tgt_mbid = mbaxid !
            mapper_Fl2a%intx_mbid = mbintxla
            mapper_Fl2a%src_context = lnd(1)%cplcompid
            mapper_Fl2a%intx_context = idintx
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fl2a%weight_identifier = wgtIdef
            mapper_Fl2a%mbname = 'mapper_Fl2a'

            if (.not. samegrid_al) then ! tri grid case
               if (iamroot_CPLID) then
                  write(logunit,*) 'iMOAB intersection between atm and land with id:', idintx
               endif
               ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mblxid, mbaxid, mbintxla)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing atm lnd intx'
                  call shr_sys_abort(subname//' ERROR in computing atm lnd intx')
               endif
#ifdef MOABDEBUG
               ! write intx only if true intx file:
               wopts = C_NULL_CHAR
               call shr_mpi_commrank( mpicom_CPLID, rank )
                  if (rank .lt. 5) then ! write only a few intx files
                  write(lnum,"(I0.2)")rank !
                  outfile = 'intx_la'//trim(lnum)// '.h5m' // C_NULL_CHAR
                  ierr = iMOAB_WriteMesh(mbintxla, outfile, wopts) ! write local intx file
                  if (ierr .ne. 0) then
                     write(logunit,*) subname,' error in writing intx file land atm '
                     call shr_sys_abort(subname//' ERROR in writing intx file ')
                  endif
               endif
#endif

               ! we also need to compute the comm graph for the second hop, from the lnd on coupler to the
               ! lnd for the intx lnd-atm context (coverage)
               !
               call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)
               type1 = 3; !  fv for lnd and atm; fv-cgll does not work anyway
               type2 = 3;

               ierr = iMOAB_ComputeCommGraph( mblxid, mbintxla, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                          lnd(1)%cplcompid, idintx)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing comm graph for second hop, ice-atm'
                  call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, ice-atm')
               endif
               ! need to compute weigths
               volumetric = 0 ! can be 1 only for FV->DGLL or FV->CGLL;

               if (atm_pg_active) then
                  dm2 = "fv"//C_NULL_CHAR
                  dofnameT="GLOBAL_ID"//C_NULL_CHAR
                  orderT = 1 !  fv-fv
               else
                  dm2 = "cgll"//C_NULL_CHAR
                  dofnameT="GLOBAL_DOFS"//C_NULL_CHAR
                  orderT = np !  it should be 4
               endif
               dm1 = "fv"//C_NULL_CHAR
               dofnameS="GLOBAL_ID"//C_NULL_CHAR
               orderS = 1  !  not much arguing
               fNoBubble = 1
               monotonicity = 0 !
               noConserve = 0
               validate = 0 ! less verbose
               fInverseDistanceMap = 0
               if (iamroot_CPLID) then
                  write(logunit,*) subname, 'launch iMOAB weights with args ', 'mbintxla=', mbintxla, ' wgtIdef=', wgtIdef, &
                     'dm1=', trim(dm1), ' orderS=',  orderS, 'dm2=', trim(dm2), ' orderT=', orderT, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT)
               endif
               ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxla, wgtIdef, &
                                                trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT) )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in iMOAB_ComputeScalarProjectionWeights lnd atm '
                  call shr_sys_abort(subname//' error in iMOAB_ComputeScalarProjectionWeights lnd atm ')
               endif

            else  ! the same mesh , atm and lnd use the same dofs, but restricted
               ! we do not compute intersection, so we will have to just send data from atm to land and viceversa, by GLOBAL_ID matching
               ! so we compute just a comm graph, between lnd and atm dofs, on the coupler; target is atm
               ! land is point cloud in this case, type1 = 2
               type1 = 3; !  full mesh for land now
               type2 = 3;  ! fv for target atm
               ierr = iMOAB_ComputeCommGraph( mblxid, mbaxid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                        lnd(1)%cplcompid, atm(1)%cplcompid)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing comm graph for second hop, lnd-atm'
                  call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, lnd-atm')
               endif
               ! context for rearrange is target in this case
               mapper_Fl2a%tgt_mbid = mbaxid
               mapper_Fl2a%intx_context = atm(1)%cplcompid

            endif ! if tri-grid
            ! we still need to define seq_flds_l2x_fields on atm cpl mesh
            if (atm_pg_active) then
               tagname = trim(seq_flds_l2x_fields)//C_NULL_CHAR
               tagtype = 1 ! dense
               numco = 1 !
               ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in defining tags for seq_flds_l2x_fields'
                  call shr_sys_abort(subname//' ERROR in coin defining tags for seq_flds_l2x_fields')
               endif
            else ! spectral case, TODO
               tagtype = 1 ! dense
            endif
         endif    ! if ((mbaxid .ge. 0) .and.  (mblxid .ge. 0) ) then
#endif
      endif ! if lnd_present
      call shr_sys_flush(logunit)

      if (lnd_c2_atm) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing mapper_Sl2a'
         endif
         call seq_map_init_rcfile(mapper_Sl2a, lnd(1), atm(1), &
            'seq_maps.rc','lnd2atm_smapname:','lnd2atm_smaptype:',samegrid_al, &
            'mapper_Sl2a initialization',esmf_map_flag)
#ifdef HAVE_MOAB
         if ((mbaxid .ge. 0) .and.  (mblxid .ge. 0) ) then
            if ( mapper_Sl2a%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Sl2a%mbname) &
                             //' mapper_Sl2a'
                endif
            endif
            mapper_Sl2a%src_mbid = mblxid
            mapper_Sl2a%tgt_mbid = mapper_Fl2a%tgt_mbid !
            mapper_Sl2a%intx_mbid = mbintxla
            mapper_Sl2a%src_context = lnd(1)%cplcompid
            mapper_Sl2a%intx_context = mapper_Fl2a%intx_context
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Sl2a%weight_identifier = wgtIdef
            mapper_Sl2a%mbname = 'mapper_Sl2a'
         endif
#endif
      endif ! if (lnd_c2_atm) then

   endif ! if atm_present

  end subroutine prep_atm_init

  !================================================================================================

  subroutine prep_atm_mrg(infodata, fractions_ax, xao_ax, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Prepare run phase, including running the merge
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_ax(:)
    type(mct_aVect)         , intent(in)    :: xao_ax(:)
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eli, eoi, eii, exi, efi, eai, emi
    type(mct_avect), pointer :: x2a_ax
    character(*), parameter  :: subname = '(prep_atm_mrg)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       ! Use fortran mod to address ensembles in merge
       eli = mod((eai-1),num_inst_lnd) + 1
       eoi = mod((eai-1),num_inst_ocn) + 1
       eii = mod((eai-1),num_inst_ice) + 1
       exi = mod((eai-1),num_inst_xao) + 1
       efi = mod((eai-1),num_inst_frc) + 1
       emi = mod((eai-1),num_inst_max) + 1

       x2a_ax => component_get_x2c_cx(atm(eai)) ! This is actually modifying x2a_ax
       call prep_atm_merge(l2x_ax(eli), o2x_ax(emi), xao_ax(exi), i2x_ax(eii), &
            fractions_ax(efi), x2a_ax)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_atm_mrg

  subroutine prep_atm_mrg_moab(infodata, xao_ax)
    use iMOAB , only : iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, &
     iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
   ! use seq_comm_mct , only : mbaxid, mbofxid ! ocean and atm-ocean flux instances
    !---------------------------------------------------------------
    ! Description
    ! Merge all ocn inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)        , pointer , intent(in)    :: xao_ax(:) ! Atm-ocn fluxes, atm grid, cpl pes; used here just for indexing

    ! Arguments
    type(mct_aVect), pointer    :: l2x_a !   needed just for indexing
    type(mct_aVect), pointer    :: o2x_a
    type(mct_aVect), pointer    :: i2x_a
    type(mct_aVect), pointer    :: xao_a
    type(mct_aVect), pointer    :: x2a_a
    ! type(mct_aVect)    :: fractions_a

    !type(mct_aVect), intent(inout) :: x2a_am
    ! ! will build x2a_am , similar to x2a_ax
    ! no averages, just one instance for atm
    ! start copy from prep_atm_merge
 !
    ! Local workspace
    real(r8) :: fracl, fraci, fraco
    integer  :: n,ka,ki,kl,ko,kx,kof,kif,klf,i,i1,o1
    integer, save :: lsize
    integer, save  :: index_x2a_Sf_lfrac, index_x2a_Sf_ifrac, index_x2a_Sf_ofrac

    character(CL),allocatable :: field_atm(:)   ! string converted to char
    character(CL),allocatable :: field_lnd(:)   ! string converted to char
    character(CL),allocatable :: field_ice(:)   ! string converted to char
    character(CL),allocatable :: field_xao(:)   ! string converted to char
    character(CL),allocatable :: field_ocn(:)   ! string converted to char
    character(CL),allocatable :: itemc_atm(:)   ! string converted to char
    character(CL),allocatable :: itemc_lnd(:)   ! string converted to char
    character(CL),allocatable :: itemc_ice(:)   ! string converted to char
    character(CL),allocatable :: itemc_xao(:)   ! string converted to char
    character(CL),allocatable :: itemc_ocn(:)   ! string converted to char
    logical :: iamroot
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    logical, save :: first_time = .true.
    type(mct_aVect_sharedindices),save :: l2x_sharedindices
    type(mct_aVect_sharedindices),save :: o2x_sharedindices
    type(mct_aVect_sharedindices),save :: i2x_sharedindices
    type(mct_aVect_sharedindices),save :: xao_sharedindices
    logical, pointer, save :: lmerge(:),imerge(:),xmerge(:),omerge(:)
    ! special for moab
    logical, pointer, save :: sharedIndex(:)
    integer, pointer, save :: lindx(:), iindx(:), oindx(:),xindx(:)
    integer, save          :: naflds, nlflds,niflds,noflds,nxflds


    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
    character(CXX) ::tagname, mct_field
    integer :: ent_type, ierr, arrsize
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
#endif

    character(*), parameter   :: subname = '(prep_atm_mrg_moab) '
    !-----------------------------------------------------------------------
    !
    call seq_comm_getdata(CPLID, iamroot=iamroot)



    if (first_time) then

      ! find out the number of local elements in moab mesh atm instance on coupler
      ! right now, we work only on FV mesh, which is a cell mesh
      ! eventually we will fix spectral case too
       ierr  = iMOAB_GetMeshInfo ( mbaxid, nvert, nvise, nbl, nsurf, nvisBC );
       if (ierr .ne. 0) then
               write(logunit,*) subname,' error in getting info '
               call shr_sys_abort(subname//' error in getting info ')
       endif
       lsize = nvise(1) ! number of active cells
       ! mct avs are used just for their fields metadata, not the actual reals
       ! (name of the fields)
       ! need these always, not only the first time
       l2x_a => l2x_ax(1)
       i2x_a => i2x_ax(1)
       o2x_a => o2x_ax(1)
       xao_a => xao_ax(1)
       x2a_a => component_get_x2c_cx(atm(1))
       naflds = mct_aVect_nRattr(x2a_a)
       nlflds = mct_aVect_nRattr(l2x_a)
       niflds = mct_aVect_nRattr(i2x_a)
       noflds = mct_aVect_nRattr(o2x_a)
       nxflds = mct_aVect_nRattr(xao_a)
       index_x2a_Sf_lfrac = mct_aVect_indexRA(x2a_a,'Sf_lfrac')
       index_x2a_Sf_ifrac = mct_aVect_indexRA(x2a_a,'Sf_ifrac')
       index_x2a_Sf_ofrac = mct_aVect_indexRA(x2a_a,'Sf_ofrac')

       !ngflds = mct_aVect_nRattr(g2x_o)
       allocate(fractions_am(lsize,5)) ! there are 5 fractions 'afrac:ifrac:ofrac:lfrac:lfrin'
       allocate(x2a_am (lsize, naflds))
       allocate(o2x_am (lsize, noflds))
       allocate(i2x_am (lsize, niflds))
       allocate(l2x_am (lsize, nlflds))
       !allocate(r2x_om (lsize, nrflds))
       allocate(xao_am (lsize, nxflds))

       allocate (sharedIndex(naflds))

       allocate(lindx(naflds), lmerge(naflds))
       allocate(iindx(naflds), imerge(naflds))
       allocate(xindx(naflds), xmerge(naflds))
       allocate(oindx(naflds), omerge(naflds))
       allocate(field_atm(naflds), itemc_atm(naflds))
       allocate(field_lnd(nlflds), itemc_lnd(nlflds))
       allocate(field_ice(niflds), itemc_ice(niflds))
       allocate(field_ocn(noflds), itemc_ocn(noflds))
       allocate(field_xao(nxflds), itemc_xao(nxflds))
       allocate(mrgstr(naflds))

       sharedIndex(:) = .false.  ! shared indices will not be set to 0 after getting them

       lindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       oindx(:) = 0
       lmerge(:)  = .true.
       imerge(:)  = .true.
       xmerge(:)  = .true.
       omerge(:)  = .true.

       do ka = 1,naflds
          field_atm(ka) = mct_aVect_getRList2c(ka, x2a_a)
          itemc_atm(ka) = trim(field_atm(ka)(scan(field_atm(ka),'_'):))
       enddo
       do kl = 1,nlflds
          field_lnd(kl) = mct_aVect_getRList2c(kl, l2x_a)
          itemc_lnd(kl) = trim(field_lnd(kl)(scan(field_lnd(kl),'_'):))
       enddo
       do ki = 1,niflds
          field_ice(ki) = mct_aVect_getRList2c(ki, i2x_a)
          itemc_ice(ki) = trim(field_ice(ki)(scan(field_ice(ki),'_'):))
       enddo
       do ko = 1,noflds
          field_ocn(ko) = mct_aVect_getRList2c(ko, o2x_a)
          itemc_ocn(ko) = trim(field_ocn(ko)(scan(field_ocn(ko),'_'):))
       enddo
       do kx = 1,nxflds
          field_xao(kx) = mct_aVect_getRList2c(kx, xao_a)
          itemc_xao(kx) = trim(field_xao(kx)(scan(field_xao(kx),'_'):))
       enddo

       call mct_aVect_setSharedIndices(l2x_a, x2a_a, l2x_SharedIndices)
       call mct_aVect_setSharedIndices(o2x_a, x2a_a, o2x_SharedIndices)
       call mct_aVect_setSharedIndices(i2x_a, x2a_a, i2x_SharedIndices)
       call mct_aVect_setSharedIndices(xao_a, x2a_a, xao_SharedIndices)

       do i=1,l2x_SharedIndices%shared_real%num_indices
          o1=l2x_SharedIndices%shared_real%aVindices2(i)
          sharedIndex(o1) = .true.
       enddo
       do i=1,o2x_SharedIndices%shared_real%num_indices
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          sharedIndex(o1) = .true.
       enddo
       do i=1,i2x_SharedIndices%shared_real%num_indices
          o1=i2x_SharedIndices%shared_real%aVindices2(i)
          sharedIndex(o1) = .true.
       enddo
       do i=1,xao_SharedIndices%shared_real%num_indices
          o1=xao_SharedIndices%shared_real%aVindices2(i)
          sharedIndex(o1) = .true.
       enddo

       ! Field naming rules
       ! Only atm states that are Sx_... will be merged
       ! Only fluxes that are F??x_... will be merged
       ! All fluxes will be multiplied by corresponding component fraction

       do ka = 1,naflds
          !--- document merge ---
          mrgstr(ka) = subname//'x2a%'//trim(field_atm(ka))//' ='
          if (field_atm(ka)(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through
          endif
          if (field_atm(ka)(1:1) == 'S' .and. field_atm(ka)(2:2) /= 'x') then
             cycle ! any state fields that are not Sx_ will just be copied
          endif

          do kl = 1,nlflds
             if (trim(itemc_atm(ka)) == trim(itemc_lnd(kl))) then
                if ((trim(field_atm(ka)) == trim(field_lnd(kl)))) then
                   if (field_lnd(kl)(1:1) == 'F') lmerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (lindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple kl field matches for ',trim(itemc_lnd(kl))
                   call shr_sys_abort(subname//' ERROR multiple kl field matches')
                endif
                lindx(ka) = kl
             endif
          end do
          do ki = 1,niflds
             if (field_ice(ki)(1:1) == 'F' .and. field_ice(ki)(2:4) == 'ioi') then
                cycle ! ignore all fluxes that are ice/ocn fluxes
             endif
             if (trim(itemc_atm(ka)) == trim(itemc_ice(ki))) then
                if ((trim(field_atm(ka)) == trim(field_ice(ki)))) then
                   if (field_ice(ki)(1:1) == 'F') imerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (iindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ki field matches for ',trim(itemc_ice(ki))
                   call shr_sys_abort(subname//' ERROR multiple ki field matches')
                endif
                iindx(ka) = ki
             endif
          end do
          do kx = 1,nxflds
             if (trim(itemc_atm(ka)) == trim(itemc_xao(kx))) then
                if ((trim(field_atm(ka)) == trim(field_xao(kx)))) then
                   if (field_xao(kx)(1:1) == 'F') xmerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (xindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple kx field matches for ',trim(itemc_xao(kx))
                   call shr_sys_abort(subname//' ERROR multiple kx field matches')
                endif
                xindx(ka) = kx
             endif
          end do
          do ko = 1,noflds
             if (trim(itemc_atm(ka)) == trim(itemc_ocn(ko))) then
                if ((trim(field_atm(ka)) == trim(field_ocn(ko)))) then
                   if (field_ocn(ko)(1:1) == 'F') omerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (oindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ko field matches for ',trim(itemc_ocn(ko))
                   call shr_sys_abort(subname//' ERROR multiple ko field matches')
                endif
                oindx(ka) = ko
             endif
          end do

          ! --- add some checks ---

          ! --- make sure all terms agree on merge or non-merge aspect ---
          if (oindx(ka) > 0 .and. xindx(ka) > 0) then
             write(logunit,*) subname,' ERROR: oindx and xindx both non-zero, not allowed ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR oindx and xindx both non-zero')
          endif

          ! --- make sure all terms agree on merge or non-merge aspect ---
          if (lindx(ka) > 0 .and. iindx(ka) > 0 .and. (lmerge(ka) .neqv. imerge(ka))) then
             write(logunit,*) subname,' ERROR: lindx and iindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR lindx and iindx merge logic error')
          endif
          if (lindx(ka) > 0 .and. xindx(ka) > 0 .and. (lmerge(ka) .neqv. xmerge(ka))) then
             write(logunit,*) subname,' ERROR: lindx and xindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR lindx and xindx merge logic error')
          endif
          if (lindx(ka) > 0 .and. oindx(ka) > 0 .and. (lmerge(ka) .neqv. omerge(ka))) then
             write(logunit,*) subname,' ERROR: lindx and oindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR lindx and oindx merge logic error')
          endif
          if (xindx(ka) > 0 .and. iindx(ka) > 0 .and. (xmerge(ka) .neqv. imerge(ka))) then
             write(logunit,*) subname,' ERROR: xindx and iindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR xindx and iindx merge logic error')
          endif
          if (xindx(ka) > 0 .and. oindx(ka) > 0 .and. (xmerge(ka) .neqv. omerge(ka))) then
             write(logunit,*) subname,' ERROR: xindx and oindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR xindx and oindx merge logic error')
          endif
          if (iindx(ka) > 0 .and. oindx(ka) > 0 .and. (imerge(ka) .neqv. omerge(ka))) then
             write(logunit,*) subname,' ERROR: iindx and oindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR iindx and oindx merge logic error')
          endif

       end do
    endif

    ! Zero attribute vector

    !call mct_avect_zero(x2a_a) ?

    !x2a_am = 0._r8
    ent_type = 1 ! cells
    tagname = trim(seq_flds_x2a_fields)//C_NULL_CHAR
    arrsize = naflds * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, x2a_am)
    if (ierr .ne. 0) then
        call shr_sys_abort(subname//' error in getting moab tags with 0 ')
    endif
    ! zero out only indices that are not shared
    do ka = 1,naflds
       if (.not. sharedIndex(ka)) then
          x2a_am(:,ka) = 0
       endif
    enddo
    ! Update surface fractions
    !    fraclist_a = 'afrac:ifrac:ofrac:lfrac:lfrin'
    kif = 2 ! kif=mct_aVect_indexRA(fractions_a,"ifrac")
    klf = 4 ! klf=mct_aVect_indexRA(fractions_a,"lfrac")
    kof = 3 ! kof=mct_aVect_indexRA(fractions_a,"ofrac")
    ! lsize = mct_avect_lsize(x2a_a)


    ! fill with fractions from atm instance

    tagname = 'afrac:ifrac:ofrac:lfrac:lfrin'//C_NULL_CHAR
    arrsize = 5 * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, fractions_am)
    if (ierr .ne. 0) then
         call shr_sys_abort(subname//' error in getting fractions_am from atm instance ')
    endif

    tagname = trim(seq_flds_o2x_fields)//C_NULL_CHAR
    arrsize = noflds * lsize !        allocate (o2x_am (lsize, noflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, o2x_am)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting o2x_am array ')
    endif

    tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
    arrsize = niflds * lsize !        allocate (i2x_am (lsize, niflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, i2x_am)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting i2x_am array ')
    endif

    tagname = trim(seq_flds_l2x_fields)//C_NULL_CHAR
    arrsize = nlflds * lsize !        allocate (l2x_am (lsize, nlflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, l2x_am)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting l2x_am array ')
    endif

    tagname = trim(seq_flds_xao_fields)//C_NULL_CHAR
    arrsize = nxflds * lsize !        allocate (xao_am (lsize, nxflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, xao_am)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting xao_om array ')
    endif


    do n = 1,lsize
       x2a_am(n, index_x2a_Sf_lfrac) = fractions_am(n, klf) ! x2a_a%rAttr(index_x2a_Sf_lfrac,n) = fractions_a%Rattr(klf,n)
       x2a_am(n, index_x2a_Sf_ifrac) = fractions_am(n, kif) ! x2a_a%rAttr(index_x2a_Sf_ifrac,n) = fractions_a%Rattr(kif,n)
       x2a_am(n, index_x2a_Sf_ofrac) = fractions_am(n, kof) ! x2a_a%rAttr(index_x2a_Sf_ofrac,n) = fractions_a%Rattr(kof,n)
    end do

    !--- document fraction operations ---
    if (first_time) then
       mrgstr(index_x2a_sf_lfrac) = trim(mrgstr(index_x2a_sf_lfrac))//' = fractions_a%lfrac'
       mrgstr(index_x2a_sf_ifrac) = trim(mrgstr(index_x2a_sf_ifrac))//' = fractions_a%ifrac'
       mrgstr(index_x2a_sf_ofrac) = trim(mrgstr(index_x2a_sf_ofrac))//' = fractions_a%ofrac'
    endif

    ! Copy attributes that do not need to be merged
    ! These are assumed to have the same name in
    ! (o2x_a and x2a_a) and in (l2x_a and x2a_a), etc.

    !--- document copy operations ---
    if (first_time) then
       !--- document merge ---
       do i=1,l2x_SharedIndices%shared_real%num_indices
          i1=l2x_SharedIndices%shared_real%aVindices1(i)
          o1=l2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = l2x%'//trim(field_lnd(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A6, I3)" )i1, ' mb-> ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,o2x_SharedIndices%shared_real%num_indices
          i1=o2x_SharedIndices%shared_real%aVindices1(i)
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field_ocn(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A6, I3)" )i1, ' mb-> ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,i2x_SharedIndices%shared_real%num_indices
          i1=i2x_SharedIndices%shared_real%aVindices1(i)
          o1=i2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = i2x%'//trim(field_ice(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A6, I3)" )i1, ' mb-> ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,xao_SharedIndices%shared_real%num_indices
          i1=xao_SharedIndices%shared_real%aVindices1(i)
          o1=xao_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = xao%'//trim(field_xao(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A6, I3)" )i1, ' mb-> ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
    endif

    !    call mct_aVect_copy(aVin=l2x_a, aVout=x2a_a, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=o2x_a, aVout=x2a_a, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=i2x_a, aVout=x2a_a, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=xao_a, aVout=x2a_a, vector=mct_usevector)
    ! we need to do something equivalent, to copy in a2x_am the tags from those shared indices
    ! call mct_aVect_copy(aVin=l2x_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=l2x_SharedIndices)
    !call mct_aVect_copy(aVin=o2x_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=o2x_SharedIndices)
    !call mct_aVect_copy(aVin=i2x_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=i2x_SharedIndices)
    !call mct_aVect_copy(aVin=xao_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=xao_SharedIndices)

    ! If flux to atm is coming only from the ocean (based on field being in o2x_a) -
    ! -- then scale by both ocean and ice fraction
    ! If flux to atm is coming only from the land or ice or coupler
    ! -- then do scale by fraction above

    do ka = 1,naflds
       !--- document merge ---
       if (first_time) then
          if (lindx(ka) > 0) then
             if (lmerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + lfrac*l2x%'//trim(field_lnd(lindx(ka)))
             else
                mrgstr(ka) = trim(mrgstr(ka))//' = lfrac*l2x%'//trim(field_lnd(lindx(ka)))
             endif
          endif
          if (iindx(ka) > 0) then
             if (imerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + ifrac*i2x%'//trim(field_ice(iindx(ka)))
             else
                mrgstr(ka) = trim(mrgstr(ka))//' = ifrac*i2x%'//trim(field_ice(iindx(ka)))
             endif
          endif
          if (xindx(ka) > 0) then
             if (xmerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + ofrac*xao%'//trim(field_xao(xindx(ka)))
             else
                mrgstr(ka) = trim(mrgstr(ka))//' = ofrac*xao%'//trim(field_xao(xindx(ka)))
             endif
          endif
          if (oindx(ka) > 0) then
             if (omerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + ofrac*o2x%'//trim(field_ocn(oindx(ka)))
             endif
             if (.not. omerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + (ifrac+ofrac)*o2x%'//trim(field_ocn(oindx(ka)))
             endif
          endif
       endif

       do n = 1,lsize
          fracl = fractions_am(n, klf) ! fractions_a%Rattr(klf,n)
          fraci = fractions_am(n, kif) ! fractions_a%Rattr(kif,n)
          fraco = fractions_am(n, kof) ! fractions_a%Rattr(kof,n)
          if (lindx(ka) > 0 .and. fracl > 0._r8) then
             if (lmerge(ka)) then
                x2a_am(n, ka) = x2a_am(n, ka) + l2x_am(n, lindx(ka)) * fracl ! x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + l2x_a%rAttr(lindx(ka),n) * fracl
             else
                x2a_am(n, ka) = l2x_am(n, lindx(ka)) * fracl ! x2a_a%rAttr(ka,n) = l2x_a%rAttr(lindx(ka),n) * fracl
             endif
          endif
          if (iindx(ka) > 0 .and. fraci > 0._r8) then
             if (imerge(ka)) then
                x2a_am(n, ka) = x2a_am(n, ka) + i2x_am(n, iindx(ka)) * fraci ! x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + i2x_a%rAttr(iindx(ka),n) * fraci
             else
                x2a_am(n, ka) = i2x_am(n, iindx(ka)) * fraci ! x2a_a%rAttr(ka,n) = i2x_a%rAttr(iindx(ka),n) * fraci
             endif
          endif
          if (xindx(ka) > 0 .and. fraco > 0._r8) then
             if (xmerge(ka)) then
                x2a_am(n, ka) = x2a_am(n, ka) + xao_am(n, xindx(ka)) * fraco !x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + xao_a%rAttr(xindx(ka),n) * fraco
             else
                x2a_am(n, ka) = xao_am(n, xindx(ka)) * fraco ! x2a_a%rAttr(ka,n) = xao_a%rAttr(xindx(ka),n) * fraco
             endif
          endif
          if (oindx(ka) > 0) then
             if (omerge(ka) .and. fraco > 0._r8) then
                x2a_am(n, ka) = x2a_am(n, ka) + o2x_am(n, oindx(ka)) * fraco ! x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco
             endif
             if (.not. omerge(ka)) then
                !--- NOTE: This IS using the ocean fields and ice fraction !! ---
                x2a_am(n, ka) = o2x_am(n, oindx(ka)) * fraci ! x2a_a%rAttr(ka,n) = o2x_a%rAttr(oindx(ka),n) * fraci
                x2a_am(n, ka) = x2a_am(n, ka) + o2x_am(n, oindx(ka)) * fraco ! x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco
             endif
          endif
       end do
    end do

! after we are done, set x2a_am to the mbaxid

    tagname = trim(seq_flds_x2a_fields)//C_NULL_CHAR
    arrsize = naflds * lsize
    ierr = iMOAB_SetDoubleTagStorage ( mbaxid, tagname, arrsize , ent_type, x2a_am)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in setting x2o_om array ')
    endif
#ifdef MOABCOMP
  !compare_mct_av_moab_tag(comp, attrVect, field, imoabApp, tag_name, ent_type, difference)
    x2a_a => component_get_x2c_cx(atm(1))
    ! loop over all fields in seq_flds_x2a_fields
    call mct_list_init(temp_list ,seq_flds_x2a_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! cell for atm, atm_pg_active
    if (iamroot) print *, subname, num_moab_exports, trim(seq_flds_x2a_fields)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call compare_mct_av_moab_tag(atm(1), x2a_a, mct_field,  mbaxid, tagname, ent_type, difference, first_time)
    enddo
    call mct_list_clean(temp_list)
#endif

#ifdef MOABDEBUG
    if (mboxid .ge. 0 ) then !  we are on coupler pes, for sure
     write(lnum,"(I0.2)")num_moab_exports
     outfile = 'AtmCplAftMm'//trim(lnum)//'.h5m'//C_NULL_CHAR
     wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
     ierr = iMOAB_WriteMesh(mbaxid, trim(outfile), trim(wopts))
    endif
#endif

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do ka = 1,naflds
             write(logunit,'(A)') trim(mrgstr(ka))
          enddo
       endif
       deallocate(mrgstr)
       deallocate(field_atm,itemc_atm)
       deallocate(field_lnd,itemc_lnd)
       deallocate(field_ice,itemc_ice)
       deallocate(field_ocn,itemc_ocn)
       deallocate(field_xao,itemc_xao)
    endif

    first_time = .false.
    ! end copy from prep_atm_merge
  end subroutine prep_atm_mrg_moab
  !================================================================================================

  subroutine prep_atm_merge( l2x_a, o2x_a, xao_a, i2x_a, fractions_a, x2a_a )

    !-----------------------------------------------------------------------
    !
    ! Arguments
    type(mct_aVect), intent(in)    :: l2x_a
    type(mct_aVect), intent(in)    :: o2x_a
    type(mct_aVect), intent(in)    :: xao_a
    type(mct_aVect), intent(in)    :: i2x_a
    type(mct_aVect), intent(in)    :: fractions_a
    type(mct_aVect), intent(inout) :: x2a_a
    !
    ! Local workspace
    real(r8) :: fracl, fraci, fraco
    integer  :: n,ka,ki,kl,ko,kx,kof,kif,klf,i,i1,o1
    integer  :: lsize
    integer  :: index_x2a_Sf_lfrac
    integer  :: index_x2a_Sf_ifrac
    integer  :: index_x2a_Sf_ofrac
    character(CL),allocatable :: field_atm(:)   ! string converted to char
    character(CL),allocatable :: field_lnd(:)   ! string converted to char
    character(CL),allocatable :: field_ice(:)   ! string converted to char
    character(CL),allocatable :: field_xao(:)   ! string converted to char
    character(CL),allocatable :: field_ocn(:)   ! string converted to char
    character(CL),allocatable :: itemc_atm(:)   ! string converted to char
    character(CL),allocatable :: itemc_lnd(:)   ! string converted to char
    character(CL),allocatable :: itemc_ice(:)   ! string converted to char
    character(CL),allocatable :: itemc_xao(:)   ! string converted to char
    character(CL),allocatable :: itemc_ocn(:)   ! string converted to char
    logical :: iamroot
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    logical, save :: first_time = .true.
    type(mct_aVect_sharedindices),save :: l2x_sharedindices
    type(mct_aVect_sharedindices),save :: o2x_sharedindices
    type(mct_aVect_sharedindices),save :: i2x_sharedindices
    type(mct_aVect_sharedindices),save :: xao_sharedindices
    logical, pointer, save :: lmerge(:),imerge(:),xmerge(:),omerge(:)
    integer, pointer, save :: lindx(:), iindx(:), oindx(:),xindx(:)
    integer, save          :: naflds, nlflds,niflds,noflds,nxflds
    character(*), parameter   :: subname = '(prep_atm_merge) '
    !-----------------------------------------------------------------------
    !
    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then

       naflds = mct_aVect_nRattr(x2a_a)
       nlflds = mct_aVect_nRattr(l2x_a)
       niflds = mct_aVect_nRattr(i2x_a)
       noflds = mct_aVect_nRattr(o2x_a)
       nxflds = mct_aVect_nRattr(xao_a)

       allocate(lindx(naflds), lmerge(naflds))
       allocate(iindx(naflds), imerge(naflds))
       allocate(xindx(naflds), xmerge(naflds))
       allocate(oindx(naflds), omerge(naflds))
       allocate(field_atm(naflds), itemc_atm(naflds))
       allocate(field_lnd(nlflds), itemc_lnd(nlflds))
       allocate(field_ice(niflds), itemc_ice(niflds))
       allocate(field_ocn(noflds), itemc_ocn(noflds))
       allocate(field_xao(nxflds), itemc_xao(nxflds))
       allocate(mrgstr(naflds))

       lindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       oindx(:) = 0
       lmerge(:)  = .true.
       imerge(:)  = .true.
       xmerge(:)  = .true.
       omerge(:)  = .true.

       do ka = 1,naflds
          field_atm(ka) = mct_aVect_getRList2c(ka, x2a_a)
          itemc_atm(ka) = trim(field_atm(ka)(scan(field_atm(ka),'_'):))
       enddo
       do kl = 1,nlflds
          field_lnd(kl) = mct_aVect_getRList2c(kl, l2x_a)
          itemc_lnd(kl) = trim(field_lnd(kl)(scan(field_lnd(kl),'_'):))
       enddo
       do ki = 1,niflds
          field_ice(ki) = mct_aVect_getRList2c(ki, i2x_a)
          itemc_ice(ki) = trim(field_ice(ki)(scan(field_ice(ki),'_'):))
       enddo
       do ko = 1,noflds
          field_ocn(ko) = mct_aVect_getRList2c(ko, o2x_a)
          itemc_ocn(ko) = trim(field_ocn(ko)(scan(field_ocn(ko),'_'):))
       enddo
       do kx = 1,nxflds
          field_xao(kx) = mct_aVect_getRList2c(kx, xao_a)
          itemc_xao(kx) = trim(field_xao(kx)(scan(field_xao(kx),'_'):))
       enddo

       call mct_aVect_setSharedIndices(l2x_a, x2a_a, l2x_SharedIndices)
       call mct_aVect_setSharedIndices(o2x_a, x2a_a, o2x_SharedIndices)
       call mct_aVect_setSharedIndices(i2x_a, x2a_a, i2x_SharedIndices)
       call mct_aVect_setSharedIndices(xao_a, x2a_a, xao_SharedIndices)

       ! Field naming rules
       ! Only atm states that are Sx_... will be merged
       ! Only fluxes that are F??x_... will be merged
       ! All fluxes will be multiplied by corresponding component fraction

       do ka = 1,naflds
          !--- document merge ---
          mrgstr(ka) = subname//'x2a%'//trim(field_atm(ka))//' ='
          if (field_atm(ka)(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through
          endif
          if (field_atm(ka)(1:1) == 'S' .and. field_atm(ka)(2:2) /= 'x') then
             cycle ! any state fields that are not Sx_ will just be copied
          endif

          do kl = 1,nlflds
             if (trim(itemc_atm(ka)) == trim(itemc_lnd(kl))) then
                if ((trim(field_atm(ka)) == trim(field_lnd(kl)))) then
                   if (field_lnd(kl)(1:1) == 'F') lmerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (lindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple kl field matches for ',trim(itemc_lnd(kl))
                   call shr_sys_abort(subname//' ERROR multiple kl field matches')
                endif
                lindx(ka) = kl
             endif
          end do
          do ki = 1,niflds
             if (field_ice(ki)(1:1) == 'F' .and. field_ice(ki)(2:4) == 'ioi') then
                cycle ! ignore all fluxes that are ice/ocn fluxes
             endif
             if (trim(itemc_atm(ka)) == trim(itemc_ice(ki))) then
                if ((trim(field_atm(ka)) == trim(field_ice(ki)))) then
                   if (field_ice(ki)(1:1) == 'F') imerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (iindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ki field matches for ',trim(itemc_ice(ki))
                   call shr_sys_abort(subname//' ERROR multiple ki field matches')
                endif
                iindx(ka) = ki
             endif
          end do
          do kx = 1,nxflds
             if (trim(itemc_atm(ka)) == trim(itemc_xao(kx))) then
                if ((trim(field_atm(ka)) == trim(field_xao(kx)))) then
                   if (field_xao(kx)(1:1) == 'F') xmerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (xindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple kx field matches for ',trim(itemc_xao(kx))
                   call shr_sys_abort(subname//' ERROR multiple kx field matches')
                endif
                xindx(ka) = kx
             endif
          end do
          do ko = 1,noflds
             if (trim(itemc_atm(ka)) == trim(itemc_ocn(ko))) then
                if ((trim(field_atm(ka)) == trim(field_ocn(ko)))) then
                   if (field_ocn(ko)(1:1) == 'F') omerge(ka) = .false.
                endif
                ! --- make sure only one field matches ---
                if (oindx(ka) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ko field matches for ',trim(itemc_ocn(ko))
                   call shr_sys_abort(subname//' ERROR multiple ko field matches')
                endif
                oindx(ka) = ko
             endif
          end do

          ! --- add some checks ---

          ! --- make sure all terms agree on merge or non-merge aspect ---
          if (oindx(ka) > 0 .and. xindx(ka) > 0) then
             write(logunit,*) subname,' ERROR: oindx and xindx both non-zero, not allowed ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR oindx and xindx both non-zero')
          endif

          ! --- make sure all terms agree on merge or non-merge aspect ---
          if (lindx(ka) > 0 .and. iindx(ka) > 0 .and. (lmerge(ka) .neqv. imerge(ka))) then
             write(logunit,*) subname,' ERROR: lindx and iindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR lindx and iindx merge logic error')
          endif
          if (lindx(ka) > 0 .and. xindx(ka) > 0 .and. (lmerge(ka) .neqv. xmerge(ka))) then
             write(logunit,*) subname,' ERROR: lindx and xindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR lindx and xindx merge logic error')
          endif
          if (lindx(ka) > 0 .and. oindx(ka) > 0 .and. (lmerge(ka) .neqv. omerge(ka))) then
             write(logunit,*) subname,' ERROR: lindx and oindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR lindx and oindx merge logic error')
          endif
          if (xindx(ka) > 0 .and. iindx(ka) > 0 .and. (xmerge(ka) .neqv. imerge(ka))) then
             write(logunit,*) subname,' ERROR: xindx and iindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR xindx and iindx merge logic error')
          endif
          if (xindx(ka) > 0 .and. oindx(ka) > 0 .and. (xmerge(ka) .neqv. omerge(ka))) then
             write(logunit,*) subname,' ERROR: xindx and oindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR xindx and oindx merge logic error')
          endif
          if (iindx(ka) > 0 .and. oindx(ka) > 0 .and. (imerge(ka) .neqv. omerge(ka))) then
             write(logunit,*) subname,' ERROR: iindx and oindx merge logic error ',trim(itemc_atm(ka))
             call shr_sys_abort(subname//' ERROR iindx and oindx merge logic error')
          endif

       end do
    endif

    ! Zero attribute vector

    call mct_avect_zero(x2a_a)

    ! Update surface fractions

    kif=mct_aVect_indexRA(fractions_a,"ifrac")
    klf=mct_aVect_indexRA(fractions_a,"lfrac")
    kof=mct_aVect_indexRA(fractions_a,"ofrac")
    lsize = mct_avect_lsize(x2a_a)

    index_x2a_Sf_lfrac = mct_aVect_indexRA(x2a_a,'Sf_lfrac')
    index_x2a_Sf_ifrac = mct_aVect_indexRA(x2a_a,'Sf_ifrac')
    index_x2a_Sf_ofrac = mct_aVect_indexRA(x2a_a,'Sf_ofrac')
    do n = 1,lsize
       x2a_a%rAttr(index_x2a_Sf_lfrac,n) = fractions_a%Rattr(klf,n)
       x2a_a%rAttr(index_x2a_Sf_ifrac,n) = fractions_a%Rattr(kif,n)
       x2a_a%rAttr(index_x2a_Sf_ofrac,n) = fractions_a%Rattr(kof,n)
    end do

    !--- document fraction operations ---
    if (first_time) then
       mrgstr(index_x2a_sf_lfrac) = trim(mrgstr(index_x2a_sf_lfrac))//' = fractions_a%lfrac'
       mrgstr(index_x2a_sf_ifrac) = trim(mrgstr(index_x2a_sf_ifrac))//' = fractions_a%ifrac'
       mrgstr(index_x2a_sf_ofrac) = trim(mrgstr(index_x2a_sf_ofrac))//' = fractions_a%ofrac'
    endif

    ! Copy attributes that do not need to be merged
    ! These are assumed to have the same name in
    ! (o2x_a and x2a_a) and in (l2x_a and x2a_a), etc.

    !--- document copy operations ---
    if (first_time) then
       !--- document merge ---
       do i=1,l2x_SharedIndices%shared_real%num_indices
          i1=l2x_SharedIndices%shared_real%aVindices1(i)
          o1=l2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = l2x%'//trim(field_lnd(i1))
       enddo
       do i=1,o2x_SharedIndices%shared_real%num_indices
          i1=o2x_SharedIndices%shared_real%aVindices1(i)
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field_ocn(i1))
       enddo
       do i=1,i2x_SharedIndices%shared_real%num_indices
          i1=i2x_SharedIndices%shared_real%aVindices1(i)
          o1=i2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = i2x%'//trim(field_ice(i1))
       enddo
       do i=1,xao_SharedIndices%shared_real%num_indices
          i1=xao_SharedIndices%shared_real%aVindices1(i)
          o1=xao_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = xao%'//trim(field_xao(i1))
       enddo
    endif

    !    call mct_aVect_copy(aVin=l2x_a, aVout=x2a_a, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=o2x_a, aVout=x2a_a, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=i2x_a, aVout=x2a_a, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=xao_a, aVout=x2a_a, vector=mct_usevector)
    call mct_aVect_copy(aVin=l2x_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=l2x_SharedIndices)
    call mct_aVect_copy(aVin=o2x_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=o2x_SharedIndices)
    call mct_aVect_copy(aVin=i2x_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=i2x_SharedIndices)
    call mct_aVect_copy(aVin=xao_a, aVout=x2a_a, vector=mct_usevector, sharedIndices=xao_SharedIndices)

    ! If flux to atm is coming only from the ocean (based on field being in o2x_a) -
    ! -- then scale by both ocean and ice fraction
    ! If flux to atm is coming only from the land or ice or coupler
    ! -- then do scale by fraction above

    do ka = 1,naflds
       !--- document merge ---
       if (first_time) then
          if (lindx(ka) > 0) then
             if (lmerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + lfrac*l2x%'//trim(field_lnd(lindx(ka)))
             else
                mrgstr(ka) = trim(mrgstr(ka))//' = lfrac*l2x%'//trim(field_lnd(lindx(ka)))
             endif
          endif
          if (iindx(ka) > 0) then
             if (imerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + ifrac*i2x%'//trim(field_ice(iindx(ka)))
             else
                mrgstr(ka) = trim(mrgstr(ka))//' = ifrac*i2x%'//trim(field_ice(iindx(ka)))
             endif
          endif
          if (xindx(ka) > 0) then
             if (xmerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + ofrac*xao%'//trim(field_xao(xindx(ka)))
             else
                mrgstr(ka) = trim(mrgstr(ka))//' = ofrac*xao%'//trim(field_xao(xindx(ka)))
             endif
          endif
          if (oindx(ka) > 0) then
             if (omerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + ofrac*o2x%'//trim(field_ocn(oindx(ka)))
             endif
             if (.not. omerge(ka)) then
                mrgstr(ka) = trim(mrgstr(ka))//' + (ifrac+ofrac)*o2x%'//trim(field_ocn(oindx(ka)))
             endif
          endif
       endif

       do n = 1,lsize
          fracl = fractions_a%Rattr(klf,n)
          fraci = fractions_a%Rattr(kif,n)
          fraco = fractions_a%Rattr(kof,n)
          if (lindx(ka) > 0 .and. fracl > 0._r8) then
             if (lmerge(ka)) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + l2x_a%rAttr(lindx(ka),n) * fracl
             else
                x2a_a%rAttr(ka,n) = l2x_a%rAttr(lindx(ka),n) * fracl
             endif
          endif
          if (iindx(ka) > 0 .and. fraci > 0._r8) then
             if (imerge(ka)) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + i2x_a%rAttr(iindx(ka),n) * fraci
             else
                x2a_a%rAttr(ka,n) = i2x_a%rAttr(iindx(ka),n) * fraci
             endif
          endif
          if (xindx(ka) > 0 .and. fraco > 0._r8) then
             if (xmerge(ka)) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + xao_a%rAttr(xindx(ka),n) * fraco
             else
                x2a_a%rAttr(ka,n) = xao_a%rAttr(xindx(ka),n) * fraco
             endif
          endif
          if (oindx(ka) > 0) then
             if (omerge(ka) .and. fraco > 0._r8) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco
             endif
             if (.not. omerge(ka)) then
                !--- NOTE: This IS using the ocean fields and ice fraction !! ---
                x2a_a%rAttr(ka,n) = o2x_a%rAttr(oindx(ka),n) * fraci
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco
             endif
          endif
       end do
    end do

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do ka = 1,naflds
             write(logunit,'(A)') trim(mrgstr(ka))
          enddo
       endif
       deallocate(mrgstr)
       deallocate(field_atm,itemc_atm)
       deallocate(field_lnd,itemc_lnd)
       deallocate(field_ice,itemc_ice)
       deallocate(field_ocn,itemc_ocn)
       deallocate(field_xao,itemc_xao)
    endif

    first_time = .false.

  end subroutine prep_atm_merge

  !================================================================================================

  subroutine prep_atm_calc_o2x_ax(fractions_ox, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_ax (note that o2x_ax is a local module variable)
#ifdef MOABDEBUG
    use iMOAB, only :  iMOAB_WriteMesh
#endif
    ! Arguments
    type(mct_aVect) , optional, intent(in) :: fractions_ox(:)
    character(len=*), optional, intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi, efi, emi
    type(mct_aVect) , pointer :: o2x_ox
    character(*), parameter   :: subname = '(prep_atm_calc_o2x_ax)'
    character(*), parameter   :: F00 = "('"//subname//" : ', 4A )"
#ifdef MOABDEBUG
    character*50             :: outfile, wopts, lnum
    integer :: ierr
#endif
    !---------------------------------------------------------------
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do emi = 1,num_inst_max
       eoi = mod((emi-1),num_inst_ocn) + 1
       efi = mod((emi-1),num_inst_frc) + 1

       o2x_ox => component_get_c2x_cx(ocn(eoi))
       if (present(fractions_ox)) then
          call seq_map_map(mapper_So2a, o2x_ox, o2x_ax(emi),&
               fldlist=seq_flds_o2x_states,norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
       else
          call seq_map_map(mapper_So2a, o2x_ox, o2x_ax(emi),&
               fldlist=seq_flds_o2x_states,norm=.true.)
       endif
       call seq_map_map(mapper_Fo2a, o2x_ox, o2x_ax(emi),&
            fldlist=seq_flds_o2x_fluxes,norm=.true.)
    enddo

    call t_drvstopf  (trim(timer))

  end subroutine prep_atm_calc_o2x_ax

  !================================================================================================

  subroutine prep_atm_calc_i2x_ax(fractions_ix, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create i2x_ax (note that i2x_ax is a local module variable)
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_ix(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eii, efi
    type(mct_aVect) , pointer :: i2x_ix
    character(*), parameter   :: subname = '(prep_atm_calc_i2x_ax)'
    character(*), parameter   :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eii = 1,num_inst_ice
       efi = mod((eii-1),num_inst_frc) + 1

       i2x_ix => component_get_c2x_cx(ice(eii))
       call seq_map_map(mapper_Si2a, i2x_ix, i2x_ax(eii), &
            fldlist=seq_flds_i2x_states, &
            avwts_s=fractions_ix(eii), avwtsfld_s='ifrac')
       call seq_map_map(mapper_Fi2a, i2x_ix, i2x_ax(eii), &
            fldlist=seq_flds_i2x_fluxes, &
            avwts_s=fractions_ix(eii), avwtsfld_s='ifrac')
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_atm_calc_i2x_ax

  !================================================================================================

  subroutine prep_atm_calc_l2x_ax(fractions_lx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_ax (note that l2x_ax is a local module variable)
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli, efi
    type(mct_avect), pointer :: l2x_lx
    character(*), parameter  :: subname = '(prep_atm_calc_l2x_ax)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       efi = mod((eli-1),num_inst_frc) + 1

       l2x_lx => component_get_c2x_cx(lnd(eli))
       call seq_map_map(mapper_Sl2a, l2x_lx, l2x_ax(eli), &
            fldlist=seq_flds_l2x_states, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
       call seq_map_map(mapper_Fl2a, l2x_lx, l2x_ax(eli), &
            fldlist=seq_flds_l2x_fluxes, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
    enddo

    call t_drvstopf  (trim(timer))

  end subroutine prep_atm_calc_l2x_ax

  !================================================================================================

  subroutine prep_atm_calc_z2x_ax(fractions_zx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create z2x_ax (note that z2x_ax is a local module variable)
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_zx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables

  end subroutine prep_atm_calc_z2x_ax

  !================================================================================================

  function prep_atm_get_l2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_l2x_ax(:)
    prep_atm_get_l2x_ax => l2x_ax(:)
  end function prep_atm_get_l2x_ax

  function prep_atm_get_i2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_i2x_ax(:)
    prep_atm_get_i2x_ax => i2x_ax(:)
  end function prep_atm_get_i2x_ax

  function prep_atm_get_o2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_o2x_ax(:)
    prep_atm_get_o2x_ax => o2x_ax(:)
  end function prep_atm_get_o2x_ax

  function prep_atm_get_z2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_z2x_ax(:)
    prep_atm_get_z2x_ax => z2x_ax(:)
  end function prep_atm_get_z2x_ax

  function prep_atm_get_mapper_So2a()
    type(seq_map), pointer :: prep_atm_get_mapper_So2a
    prep_atm_get_mapper_So2a => mapper_So2a
  end function prep_atm_get_mapper_So2a

  function prep_atm_get_mapper_Fo2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fo2a
    prep_atm_get_mapper_Fo2a => mapper_Fo2a
  end function prep_atm_get_mapper_Fo2a

  function prep_atm_get_mapper_Sof2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Sof2a
    prep_atm_get_mapper_Sof2a => mapper_Sof2a
  end function prep_atm_get_mapper_Sof2a

  function prep_atm_get_mapper_Fof2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fof2a
    prep_atm_get_mapper_Fof2a => mapper_Fof2a
  end function prep_atm_get_mapper_Fof2a

  function prep_atm_get_mapper_Sl2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Sl2a
    prep_atm_get_mapper_Sl2a => mapper_Sl2a
  end function prep_atm_get_mapper_Sl2a

  function prep_atm_get_mapper_Fl2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fl2a
    prep_atm_get_mapper_Fl2a => mapper_Fl2a
  end function prep_atm_get_mapper_Fl2a

  function prep_atm_get_mapper_Si2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Si2a
    prep_atm_get_mapper_Si2a => mapper_Si2a
  end function prep_atm_get_mapper_Si2a

  function prep_atm_get_mapper_Fi2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fi2a
    prep_atm_get_mapper_Fi2a => mapper_Fi2a
  end function prep_atm_get_mapper_Fi2a

  !================================================================================================

end module prep_atm_mod
