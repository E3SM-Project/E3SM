module prep_ocn_mod

  use shr_kind_mod,     only: R8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: CX => shr_kind_CX, CXX => shr_kind_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_atm, num_inst_rof, num_inst_ice
  use seq_comm_mct,     only: num_inst_glc, num_inst_wav, num_inst_ocn
  use seq_comm_mct,     only: num_inst_xao, num_inst_frc
  use seq_comm_mct,     only: num_inst_max
  use seq_comm_mct,     only: CPLID, OCNID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs

  use seq_comm_mct,     only: mboxid ! iMOAB id for mpas ocean migrated mesh to coupler pes
  use seq_comm_mct,     only: mbrmapro ! iMOAB id for map read from rof2ocn map file
  use seq_comm_mct,     only: mbrxoid ! iMOAB id for rof on coupler in ocean context;
  use seq_comm_mct,     only: mbrxid   ! iMOAB id of moab rof read on couple pes
  use seq_comm_mct,     only : atm_pg_active  ! whether the atm uses FV mesh or not ; made true if fv_nphys > 0
  use seq_comm_mct,     only : mbaxid   ! iMOAB id for atm migrated mesh to coupler pes
  use seq_comm_mct,     only : mbintxao ! iMOAB id for intx mesh between atm and ocean
  use seq_comm_mct,     only : mbintxoa ! iMOAB id for intx mesh between ocean and atmosphere
  use seq_comm_mct,     only : mhid     ! iMOAB id for atm instance
  use seq_comm_mct,     only : mhpgid   ! iMOAB id for atm pgx grid, on atm pes; created with se and gll grids
  ! use dimensions_mod,   only : np     ! for atmosphere degree
  use seq_comm_mct,     only : mbixid   ! iMOAB for sea-ice migrated to coupler
  use seq_comm_mct,     only : CPLALLICEID
  use seq_comm_mct,     only : seq_comm_iamin
  use seq_comm_mct,     only : num_moab_exports

  use seq_comm_mct,     only : seq_comm_getinfo => seq_comm_setptrs

  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_map_type_mod
  use seq_map_mod        !  will have also moab_map_init_rcfile , seq_map_set_type
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
#ifdef MOABCOMP
  use component_type_mod, only:  compare_mct_av_moab_tag
#endif
  use component_type_mod, only: ocn, atm, ice, rof, wav, glc
  use iso_c_binding

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_ocn_init
  public :: prep_ocn_mrg
  ! moab version
  public :: prep_ocn_mrg_moab

  public :: prep_ocn_accum
  public :: prep_ocn_accum_moab
  public :: prep_ocn_accum_avg
  public :: prep_ocn_accum_avg_moab

  public :: prep_ocn_calc_a2x_ox

  public :: prep_ocn_calc_i2x_ox
  public :: prep_ocn_calc_r2x_ox
  public :: prep_ocn_calc_g2x_ox
  public :: prep_ocn_shelf_calc_g2x_ox
  public :: prep_ocn_calc_w2x_ox

  public :: prep_ocn_get_a2x_ox
  public :: prep_ocn_get_r2x_ox
  public :: prep_ocn_get_i2x_ox
  public :: prep_ocn_get_g2x_ox
  public :: prep_ocn_get_w2x_ox

  public :: prep_ocn_get_x2oacc_ox
  public :: prep_ocn_get_x2oacc_ox_cnt

  public :: prep_ocn_get_x2oacc_om ! will return a pointer to the local private matrix
  public :: prep_ocn_get_x2oacc_om_cnt
#ifdef SUMMITDEV_PGI
  ! Sarat: Dummy variable added to workaround PGI compiler bug (PGI 17.9) as of Oct 23, 2017
  public :: dummy_pgibugfix
#endif

  public :: prep_ocn_get_mapper_Sa2o
  public :: prep_ocn_get_mapper_Va2o
  public :: prep_ocn_get_mapper_Fa2o
  public :: prep_ocn_get_mapper_Fr2o
  public :: prep_ocn_get_mapper_Rr2o_liq
  public :: prep_ocn_get_mapper_Rr2o_ice
  public :: prep_ocn_get_mapper_SFi2o
  public :: prep_ocn_get_mapper_Rg2o_liq
  public :: prep_ocn_get_mapper_Rg2o_ice
  public :: prep_ocn_get_mapper_Sg2o
  public :: prep_ocn_get_mapper_Fg2o
  public :: prep_ocn_get_mapper_Sw2o

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_ocn_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sa2o
  type(seq_map), pointer :: mapper_Va2o
  type(seq_map), pointer :: mapper_Fa2o
  type(seq_map), pointer :: mapper_Fr2o
  type(seq_map), pointer :: mapper_Rr2o_liq
  type(seq_map), pointer :: mapper_Rr2o_ice
  type(seq_map), pointer :: mapper_SFi2o
  type(seq_map), pointer :: mapper_Rg2o_liq
  type(seq_map), pointer :: mapper_Rg2o_ice
  type(seq_map), pointer :: mapper_Fg2o
  type(seq_map), pointer :: mapper_Sg2o
  type(seq_map), pointer :: mapper_Sw2o

  ! attribute vectors
  type(mct_aVect), pointer :: a2x_ox(:) ! Atm export, ocn grid, cpl pes
  type(mct_aVect), pointer :: r2x_ox(:) ! Rof export, ocn grid, cpl pes
  type(mct_aVect), pointer :: i2x_ox(:) ! Ice export, ocn grid, cpl pes
  type(mct_aVect), pointer :: g2x_ox(:) ! Glc export, ocn grid, cpl pes
  type(mct_aVect), pointer :: w2x_ox(:) ! Wav export, ocn grid, cpl pes

  type(mct_aVect), target  :: x2o_ox_inst  ! multi instance for averaging

  ! accumulation variables
  type(mct_aVect), pointer :: x2oacc_ox(:)  ! Ocn import, ocn grid, cpl pes
  integer        , target  :: x2oacc_ox_cnt ! x2oacc_ox: number of time samples accumulated

  ! accumulation variables for moab data
  real (kind=R8) , allocatable, private, target :: x2oacc_om (:,:)   ! Ocn import, ocn grid, cpl pes, moab array
  integer        , target  :: x2oacc_om_cnt ! x2oacc_ox: number of time samples accumulated, in moab array
  integer                  :: arrSize_x2o_om !   this will be a module variable, size moabLocal_size * nof

  ! other module variables
  integer       :: mpicom_CPLID   ! MPI cpl communicator
  logical       :: iamroot_CPLID  ! .true. => CPLID masterproc
  logical       :: flood_present  ! .true.  => rof is computing flood
  character(CS) :: vect_map       ! vector mapping type
  logical       :: x2o_average    ! logical for x2o averaging to 1 ocean instance from multi instances
#ifdef SUMMITDEV_PGI
  ! Sarat: Dummy variable added to workaround PGI compiler bug (PGI 17.9) as of Oct 23, 2017
  logical       :: dummy_pgibugfix
#endif
  !================================================================================================


  real (kind=R8) , allocatable, private :: fractions_om (:,:) ! will retrieve the fractions from ocean, and use them
  !  they were init with
  ! character(*),parameter :: fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad' in moab, on the fractions
  real (kind=R8) , allocatable, private :: x2o_om (:,:)
  real (kind=R8) , allocatable, private :: a2x_om (:,:)
  real (kind=R8) , allocatable, private :: i2x_om (:,:)
  real (kind=R8) , allocatable, private :: r2x_om (:,:)
  real (kind=R8) , allocatable, private :: xao_om (:,:)

  ! this will be constructed first time, and be used to copy fields for shared indices
  ! between xao and x2o
  character(CXX) :: shared_fields_xao_x2o
  ! will need some array to hold the data for copying
  real(R8) , allocatable, save :: shared_values(:) ! will be the size of shared indices * lsize
  integer    :: size_of_shared_values

  logical                  :: iamin_CPLALLICEID     ! pe associated with CPLALLICEID
contains

  !================================================================================================


   subroutine print_weight_map_details(subname, mbintxao, maptype, identifier, &
      srcmethod, srcorder, srcdofname, tgtmethod, tgtorder, tgtdofname, fvmethod, nobubble, &
      monotonicity, volumetric, inversemap, noconserve, validate)

    integer, intent(in)      :: mbintxao, nobubble, volumetric, inversemap, noconserve, validate
    integer, intent(in)      :: srcorder, tgtorder, monotonicity
    character(*), intent(in) :: subname, maptype, identifier, srcmethod, tgtmethod, &
                                srcdofname, tgtdofname, fvmethod

     write(logunit,*) subname, ': iMOAB computing remapping weights', maptype, &
            ' and identifier ', identifier, ' with arguments: mbintxao=', mbintxao,  &
            ', source (method/order/doftag)=', trim(srcmethod), srcorder, trim(srcdofname), &
            ', target (method/order/doftag)=', trim(tgtmethod), tgtorder, trim(tgtdofname), &
            ', fvMethod=', trim(fvmethod), ', nobubble=', nobubble, &
            ', monotonicity=', monotonicity, ', volumetric=', volumetric, &
            ', fInverseDistanceMap=', inversemap, ', noConserve=', noconserve, &
            ', validate=', validate

   end subroutine print_weight_map_details


  subroutine prep_ocn_init(infodata, atm_c2_ocn, atm_c2_ice, ice_c2_ocn, rof_c2_ocn, &
       wav_c2_ocn, glc_c2_ocn, glcshelf_c2_ocn)

    use iMOAB, only: iMOAB_ComputeMeshIntersectionOnSphere, iMOAB_RegisterApplication, &
      iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_ComputeCommGraph, iMOAB_ComputeScalarProjectionWeights, &
      iMOAB_MigrateMapMesh, iMOAB_WriteLocalMesh, iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorage, &
      iMOAB_WriteMappingWeightsToFile, iMOAB_SetMapGhostLayers
    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables except for accumulators
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_ocn ! .true.=>atm to ocn coupling on
    logical                 , intent(in)    :: atm_c2_ice ! .true.=>atm to ice coupling on
    logical                 , intent(in)    :: ice_c2_ocn ! .true.=>ice to ocn coupling on
    logical                 , intent(in)    :: rof_c2_ocn ! .true.=>rof to ocn coupling on
    logical                 , intent(in)    :: wav_c2_ocn ! .true.=>wav to ocn coupling on
    logical                 , intent(in)    :: glc_c2_ocn ! .true.=>glc to ocn coupling on
    logical                 , intent(in)    :: glcshelf_c2_ocn ! .true.=>glc ice shelf to ocn coupling on
    !
    ! Local Variables
    logical                  :: esmf_map_flag  ! .true. => use esmf for mapping
    logical                  :: ocn_present    ! .true.  => ocn is present
    logical                  :: atm_present    ! .true.  => atm is present
    logical                  :: ice_present    ! .true.  => ice is present
    logical                  :: samegrid_ao    ! samegrid atm and ocean
    logical                  :: samegrid_og    ! samegrid glc and ocean
    logical                  :: samegrid_ow    ! samegrid ocean and wave
    logical                  :: samegrid_ro    ! samegrid runoff and ocean
    integer                  :: atm_nx, atm_ny
    integer                  :: lsize_o
    integer                  :: egi, eri
    integer                  :: ewi, eai, eii, eoi
    character(CL)            :: ocn_gnam       ! ocn grid
    character(CL)            :: atm_gnam       ! atm grid
    character(CL)            :: rof_gnam       ! rof grid
    character(CL)            :: wav_gnam       ! wav grid
    character(CL)            :: glc_gnam       ! glc grid
    type(mct_avect), pointer :: o2x_ox
    type(mct_avect), pointer :: x2o_ox
    character(*), parameter  :: subname = '(prep_ocn_init)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    character(*), parameter  :: F01 = "('"//subname//" : ', A, I8 )"

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

    integer                  :: rmapid, rmapid2  ! external id to identify the moab app ; 2 is for rof in ocean context (coverage)
    integer                  :: type_grid !
    integer                  :: context_id, direction
    character*32             :: prefix_output ! for writing a coverage file for debugging
    integer                  :: rank_on_cpl ! just for debugging
! these are just to zero out r2x fields on ocean
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
    integer  mlsize ! moab local ocean size
    integer  nrflds  ! number of rof fields projected on land
    integer arrsize  ! for setting the r2x fields on land to 0
    integer ent_type ! for setting tags
    integer noflds   ! used for number of fields in allocating moab accumulated array x2oacc_om
    real (kind=R8) , allocatable :: tmparray (:) ! used to set the r2x fields to 0
    integer  nghlay ! used to set the number of ghost layers, needed for bilinear map
    integer  nghlay_tgt

    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         ocn_present=ocn_present       , &
         atm_present=atm_present       , &
         ice_present=ice_present       , &
         flood_present=flood_present   , &
         vect_map=vect_map             , &
         atm_gnam=atm_gnam             , &
         ocn_gnam=ocn_gnam             , &
         rof_gnam=rof_gnam             , &
         wav_gnam=wav_gnam             , &
         atm_nx=atm_nx                 , &
         atm_ny=atm_ny                 , &
         glc_gnam=glc_gnam             , &
         esmf_map_flag=esmf_map_flag   )

    allocate(mapper_Sa2o)
    allocate(mapper_Va2o)
    allocate(mapper_Fa2o)
    allocate(mapper_Fr2o)
    allocate(mapper_Rr2o_liq)
    allocate(mapper_Rr2o_ice)
    allocate(mapper_SFi2o)
    allocate(mapper_Rg2o_liq)
    allocate(mapper_Rg2o_ice)
    allocate(mapper_Sg2o)
    allocate(mapper_Fg2o)
    allocate(mapper_Sw2o)

    if (ocn_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       o2x_ox => component_get_c2x_cx(ocn(1))
       x2o_ox => component_get_x2c_cx(ocn(1))
       lsize_o = mct_aVect_lsize(o2x_ox)

       ! x2o_average setup logic
       if (num_inst_max == num_inst_ocn) then
          ! standard multi-instance merge
          x2o_average = .false.
       elseif (num_inst_max > 1 .and. num_inst_ocn == 1) then
          ! averaging ocean merge
          x2o_average = .true.
          if (iamroot_CPLID) then
             write(logunit,F01) 'x2o averaging on over instances =',num_inst_max
          end if
          call mct_aVect_init(x2o_ox_inst, x2o_ox, lsize_o)
          call mct_aVect_zero(x2o_ox_inst)
       else
          ! not allowed
          write(logunit,F00) ' ERROR in x2o_average setup logic '
          call shr_sys_abort(subname//' ERROR in x2o_average setup logic')
       endif

       allocate(a2x_ox(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_ox(eai), rList=seq_flds_a2x_fields, lsize=lsize_o)
          call mct_aVect_zero(a2x_ox(eai))
       enddo
       allocate(r2x_ox(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_aVect_init(r2x_ox(eri), rList=seq_flds_r2x_fields, lsize=lsize_o)
          call mct_aVect_zero(r2x_ox(eri))
       enddo
       allocate(g2x_ox(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(g2x_ox(egi), rList=seq_flds_g2x_fields, lsize=lsize_o)
          call mct_aVect_zero(g2x_ox(egi))
       end do
       allocate(w2x_ox(num_inst_wav))
       do ewi = 1,num_inst_wav
          call mct_aVect_init(w2x_ox(ewi), rList=seq_flds_w2x_fields, lsize=lsize_o)
          call mct_aVect_zero(w2x_ox(ewi))
       enddo
       allocate(i2x_ox(num_inst_ice))
       do eii = 1,num_inst_ice
          call mct_aVect_init(i2x_ox(eii), rList=seq_flds_i2x_fields, lsize=lsize_o)
          call mct_aVect_zero(i2x_ox(eii))
       enddo

       allocate(x2oacc_ox(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_avect_init(x2oacc_ox(eoi), x2o_ox, lsize_o)
          call mct_aVect_zero(x2oacc_ox(eoi))
       end do
       x2oacc_ox_cnt = 0

       ! moab accumulation variable has to be allocated here too, because restart needs it
       ! it resulted in unexpected crashes, because we were allocating it only
       ! during "first_time" entering merge routine; this was wrong
       ! allocate accumulation variable , parallel to x2o_om
       noflds = mct_aVect_nRattr(x2o_ox) ! these are saved after first time
       ! size of the x2oacc_om depends on the size of the ocean mesh locally

       ! find out the number of local elements in moab mesh ocean instance on coupler
       ierr  = iMOAB_GetMeshInfo ( mboxid, nvert, nvise, nbl, nsurf, nvisBC )
       if (ierr .ne. 0) then
          write(logunit,*) subname,' cant get size of ocn mesh'
          call shr_sys_abort(subname//' ERROR in getting size of ocn mesh')
       endif
         ! ocn is cell mesh on coupler side
       mlsize = nvise(1)
       allocate(x2oacc_om(mlsize, noflds))
       x2oacc_om_cnt = 0
       x2oacc_om(:,:)=0.

       ! moab accumulation variable

       samegrid_ao = .true.
       samegrid_ro = .true.
       samegrid_ow = .true.
       samegrid_og = .true.
       if (trim(atm_gnam) /= trim(ocn_gnam)) samegrid_ao = .false.
       if (trim(rof_gnam) /= trim(ocn_gnam)) samegrid_ro = .false.
       if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = .false.
       if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_og = .false.

       if (atm_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fa2o'
          end if
          call seq_map_init_rcfile(mapper_Fa2o, atm(1), ocn(1), &
               'seq_maps.rc','atm2ocn_fmapname:','atm2ocn_fmaptype:',samegrid_ao, &
               'mapper_Fa2o initialization',esmf_map_flag)
          call shr_sys_flush(logunit)
#ifdef HAVE_MOAB
          ! Call moab intx only if atm and ocn are init in moab
          if ((mbaxid .ge. 0) .and.  (mboxid .ge. 0)) then
            appname = "ATM_OCN_COU"//C_NULL_CHAR
            ! idintx is a unique number of MOAB app that takes care of intx between ocn and atm mesh
            idintx = 100*atm(1)%cplcompid + ocn(1)%cplcompid ! something different, to differentiate it
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, idintx, mbintxao)
            if (ierr .ne. 0) then
              write(logunit,*) subname,' error in registering atm ocn intx'
              call shr_sys_abort(subname//' ERROR in registering atm ocn intx')
            endif

            mapper_Fa2o%src_mbid = mbaxid
            mapper_Fa2o%tgt_mbid = mboxid
            mapper_Fa2o%intx_mbid = mbintxao
            mapper_Fa2o%src_context = atm(1)%cplcompid
            wgtIdef = 'scalar'//C_NULL_CHAR
            mapper_Fa2o%weight_identifier = wgtIdef
            mapper_Fa2o%mbname = 'mapper_Fa2o'

            ! we also need to compute the comm graph for the second hop, from the atm on coupler to the
            ! atm for the intx atm-ocn context (coverage)
            call seq_comm_getinfo(CPLID, mpigrp=mpigrp_CPLID)

            ! next, let us compute the ATM and OCN data transfer
            if (.not. samegrid_ao) then ! not a data OCN model

               ! for bilinear maps, we need to have a layer of ghosts on source
               nghlay = 1  ! number of ghost layers
               nghlay_tgt = 0
               ierr   = iMOAB_SetMapGhostLayers( mbintxao, nghlay, nghlay_tgt )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in setting the number of layers'
                  call shr_sys_abort(subname//' error in setting the number of layers')
               endif
               ! first compute the overlap mesh between mbaxid (ATM) and mboxid (OCN) on coupler PEs
               ierr =  iMOAB_ComputeMeshIntersectionOnSphere (mbaxid, mboxid, mbintxao)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing ATM-OCN intersection'
                  call shr_sys_abort(subname//' ERROR in computing ATM-OCN intersection')
               endif
               if (iamroot_CPLID) then
                  write(logunit,*) 'iMOAB intersection completed between atm and ocean with id:', idintx
               end if

               if (atm_pg_active) then
                  type1 = 3; ! FV for ATM; CGLL does not work correctly in parallel at the moment
               else
                  type1 = 1 ! This projection works (CGLL to FV), but reverse does not (FV - CGLL)
               endif
               type2 = 3;  ! FV mesh on coupler OCN
               ! ierr      = iMOAB_ComputeCommGraph( mboxid, mbintxoa, &mpicom_CPLID, &mpigrp_CPLID, &mpigrp_CPLID, &type1, &type2,
               !                              &ocn_id, &idintx)
               ierr = iMOAB_ComputeCommGraph( mbaxid, mbintxao, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                          atm(1)%cplcompid, idintx)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing comm graph for second hop, atm-ocn'
                  call shr_sys_abort(subname//' ERROR in computing comm graph for second hop, atm-ocn')
               endif
               ! now take care of the mapper
               if ( mapper_Fa2o%src_mbid .gt. -1 ) then
                  if (iamroot_CPLID) then
                        write(logunit,F00) 'overwriting '//trim(mapper_Fa2o%mbname) &
                              //' mapper_Fa2o'
                  endif
               endif

               ! To project fields from ATM to OCN grid, we need to define
               ! ATM a2x fields to OCN grid on coupler side
               tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
               tagtype = 1 ! dense
               numco = 1 !
               ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in defining tags for seq_flds_a2x_fields on ocn cpl'
                  call shr_sys_abort(subname//' ERROR in coin defining tags for seq_flds_a2x_fields on ocn cpl')
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

               ! First compute the non-conservative bilinear map for projection of scalar fields
               if (iamroot_CPLID) then
                  call print_weight_map_details(subname, mbintxao, "FV-FV", "bilinear", &
                     trim(dm1), orderS, trim(dofnameS), trim(dm2), orderT, trim(dofnameT), "bilinear", &
                     fNoBubble, monotonicity, volumetric, fInverseDistanceMap, noConserve, validate)
               endif
               ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxao, 'bilinear'//C_NULL_CHAR, &
                                                trim(dm1), orderS, trim(dm2), orderT, 'bilin'//C_NULL_CHAR, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT) )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing ao weights '
                  call shr_sys_abort(subname//' ERROR in computing ao weights ')
               endif

               ! ierr = iMOAB_WriteMappingWeightsToFile(mbintxao, 'bilinear'//C_NULL_CHAR, 'bilinear_a2o.nc'//C_NULL_CHAR)

               ! Next compute the conservative map for projection of flux fields
               if (iamroot_CPLID) then
                  call print_weight_map_details(subname, mbintxao, "FV-FV", wgtIdef, &
                     trim(dm1), orderS, trim(dofnameS), trim(dm2), orderT, trim(dofnameT), "", &
                     fNoBubble, monotonicity, volumetric, fInverseDistanceMap, noConserve, validate)
               endif
               ierr = iMOAB_ComputeScalarProjectionWeights ( mbintxao, wgtIdef, &
                                                trim(dm1), orderS, trim(dm2), orderT, ''//C_NULL_CHAR, &
                                                fNoBubble, monotonicity, volumetric, fInverseDistanceMap, &
                                                noConserve, validate, &
                                                trim(dofnameS), trim(dofnameT) )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing ao weights '
                  call shr_sys_abort(subname//' ERROR in computing ao weights ')
               endif

               mapper_Fa2o%intx_context = idintx
               !! All done for mapper_Fa2o --

#ifdef MOABDEBUG
               wopts = C_NULL_CHAR
               call shr_mpi_commrank( mpicom_CPLID, rank )
               if (rank .lt. 5) then
                  write(lnum,"(I0.2)")rank !
                  outfile = 'intx_ao_'//trim(lnum)// '.h5m' // C_NULL_CHAR
                  ierr = iMOAB_WriteMesh(mbintxao, outfile, wopts) ! write local intx file
                  if (ierr .ne. 0) then
                     write(logunit,*) subname,' error in writing intx file '
                     call shr_sys_abort(subname//' ERROR in writing intx file ')
                  endif
               endif
#endif
            else ! if (samegrid_ao)

               ! ATM and OCN components use the same mesh and DoF numbering (OCN is a subset of ATM);
               ! We do not need to compute intersection since the "map" from ATM to OCN is essentially a
               ! permutation operator and can be achieved by simply sending/receiving data from ATM to OCN
               ! and viceversa, based on element GLOBAL_ID matching. In order to seamless produce the
               ! permutation operator, we will compute a communication graph between ATM and OCN DoFs on the
               ! coupler.
              if (atm_pg_active) then
                  type1 = 3; ! FV for ATM; CGLL does not work correctly in parallel at the moment
              else
                  type1 = 1 ! This projection works (CGLL to FV), but reverse does not (FV - CGLL)
              endif
              type2 = 3;  ! FV mesh on coupler OCN
              ierr = iMOAB_ComputeCommGraph( mbaxid, mboxid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, type1, type2, &
                                      atm(1)%cplcompid, ocn(1)%cplcompid )
              if (ierr .ne. 0) then
                write(logunit,*) subname,' error in computing communication graph for second hop, ATM-OCN'
                call shr_sys_abort(subname//' ERROR in computing communication graph for second hop, ATM-OCN')
              endif
              mapper_Fa2o%intx_context = ocn(1)%cplcompid

            endif ! if (.not. samegrid_ao)

         endif ! if ((mbaxid .ge. 0) .and. (mboxid .ge. 0))
! endif HAVE_MOAB
#endif
       end if ! if (atm_present)

       ! atm_c2_ice flag is here because ICE and OCN are constrained to be on the same
       ! grid so the ATM->ICE mapping is set to the atm->ocn mapping to improve performance
       if (atm_c2_ocn .or. atm_c2_ice) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2o'
          end if
          call seq_map_init_rcfile(mapper_Sa2o, atm(1), ocn(1), &
               'seq_maps.rc','atm2ocn_smapname:','atm2ocn_smaptype:',samegrid_ao, &
               'mapper_Sa2o initialization',esmf_map_flag)

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Va2o'
          end if
          call seq_map_init_rcfile(mapper_Va2o, atm(1), ocn(1), &
               'seq_maps.rc','atm2ocn_vmapname:','atm2ocn_vmaptype:',samegrid_ao, &
               'mapper_Va2o initialization',esmf_map_flag)

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Va2o vect with vect_map = ',trim(vect_map)
          end if
          call seq_map_initvect(mapper_Va2o, vect_map, atm(1), ocn(1), string='mapper_Va2o initvect')

          ! will use the same map for mapper_Sa2o and Va2o, using the bilinear map option
          if ((mbaxid .ge. 0) .and.  (mboxid .ge. 0)) then

            ! now take care of the 2 new mappers
            if ( mapper_Sa2o%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Sa2o%mbname) &
                             //' mapper_Sa2o'
                endif
            endif
            mapper_Sa2o%src_mbid = mbaxid
            mapper_Sa2o%tgt_mbid = mboxid
            mapper_Sa2o%intx_mbid = mbintxao
            mapper_Sa2o%src_context = atm(1)%cplcompid
            mapper_Sa2o%intx_context = mapper_Fa2o%intx_context
            mapper_Sa2o%weight_identifier = 'bilinear'//C_NULL_CHAR
            mapper_Sa2o%mbname = 'mapper_Sa2o'

            if ( mapper_Va2o%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Va2o%mbname) &
                             //' mapper_Va2o'
                endif
            endif
            mapper_Va2o%src_mbid = mbaxid
            mapper_Va2o%tgt_mbid = mboxid
            mapper_Va2o%intx_mbid = mbintxao
            mapper_Va2o%src_context = atm(1)%cplcompid
            mapper_Va2o%intx_context = mapper_Fa2o%intx_context
            mapper_Va2o%weight_identifier = 'bilinear'//C_NULL_CHAR
            mapper_Va2o%mbname = 'mapper_Va2o'

          endif ! if ((mbaxid .ge. 0) .and.  (mboxid .ge. 0))
       endif ! if (atm_c2_ocn .or. atm_c2_ice)
       call shr_sys_flush(logunit)

       ! needed for domain checking
       if (ice_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_SFi2o'
          end if
          call seq_map_init_rearrolap(mapper_SFi2o, ice(1), ocn(1), 'mapper_SFi2o')
#ifdef HAVE_MOAB
          if ( (mbixid .ge. 0) .and. (mboxid .ge. 0)) then
            ! moab also will do just a rearrange, hopefully, in this case, based on the comm graph
            !   that is computed here
            call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)   !  second group, the coupler group CPLID is global variable

            type1 = 3
            type2 = 3 ! FV-FV graph

            ! iMOAB: compute the communication graph for ICE-OCN, based on the same global id
            ! it will be a simple permutation from ice mesh directly to ocean, using the comm graph computed here
            ierr = iMOAB_ComputeCommGraph( mbixid, mboxid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, &
               type1, type2, ice(1)%cplcompid, ocn(1)%cplcompid)
            if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing graph ice - ocn x '
                  call shr_sys_abort(subname//' ERROR  in computing graph ice - ocn x ')
            endif

            ! define tags according to the seq_flds_i2x_fields
            tagtype = 1  ! dense, double
            numco = 1 !  one value per cell / entity
            tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
            if ( ierr == 1 ) then
               call shr_sys_abort( subname//' ERROR: cannot define tags for ice proj to ocn' )
            end if
            if ( mapper_SFi2o%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_SFi2o%mbname) &
                             //' mapper_SFi2o'
                endif
            endif
            mapper_SFi2o%src_mbid = mbixid
            mapper_SFi2o%tgt_mbid = mboxid
            ! no intersection, so will have to do without it
            mapper_SFi2o%src_context = ice(1)%cplcompid
            mapper_SFi2o%intx_context = ocn(1)%cplcompid
            mapper_SFi2o%mbname = 'mapper_SFi2o'

            if(mapper_SFi2o%copy_only) then
               call seq_map_set_type(mapper_SFi2o, mbixid, 1) ! type is cells
            endif

         endif

#endif

       endif ! if (ice_present)
       call shr_sys_flush(logunit)

       if (rof_c2_ocn) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rr2o_liq'
          end if
          call seq_map_init_rcfile(mapper_Rr2o_liq, rof(1), ocn(1), &
               'seq_maps.rc', 'rof2ocn_liq_rmapname:', 'rof2ocn_liq_rmaptype:',samegrid_ro, &
               'mapper_Rr2o_liq  initialization',esmf_map_flag)

#ifdef HAVE_MOAB
          appname = "ROF_OCN_COU"//CHAR(0)
            ! rmapid  is a unique external number of MOAB app that takes care of map between rof and ocn mesh
          rmapid = 100*rof(1)%cplcompid + ocn(1)%cplcompid ! something different, to differentiate it
          ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, rmapid, mbrmapro)
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in registering rof 2 ocn moab map '
             call shr_sys_abort(subname//' ERROR in registering  rof 2 ocn moab map ')
          endif
 ! integer, public :: mboxid   ! iMOAB id for mpas ocean already migrated mesh to coupler pes
          type_grid = 3 ! this is type of grid, maybe should be saved on imoab app ?
          call moab_map_init_rcfile(mbrmapro, mboxid, type_grid, rof(1), ocn(1), &
               'seq_maps.rc', 'rof2ocn_liq_rmapname:', 'rof2ocn_liq_rmaptype:',samegrid_ro, &
               'mapper_Rr2o_liq moab initialization',esmf_map_flag)
          ! this is a special rof mesh redistribution, for the ocean context
          ! it will be used to project from rof to ocean
          ! the mesh will be migrated, to be able to do the second hop
          appname = "ROF_OCOU"//C_NULL_CHAR
          ! rmapid  is a unique external number of MOAB app that identifies runoff on coupler side
          rmapid2 = 100*rof(1)%cplcompid ! this is a special case, because we also have a regular coupler instance mbrxid
          ierr = iMOAB_RegisterApplication(trim(appname), mpicom_CPLID, rmapid2, mbrxoid)
          if (ierr .ne. 0) then
             write(logunit,*) subname,' error in registering rof on coupler in ocean context '
             call shr_sys_abort(subname//' ERROR in registering  rof on coupler in ocean context ')
          endif
          ! this code was moved from prep_rof_ocn_moab, because we will do everything on coupler side, not
          ! needed to be on joint comm anymore for the second hop

      !  it read on the coupler side, from file, the scrip mosart, that has a full mesh;
      !  also migrate rof mesh on coupler pes, in ocean context, mbrxoid (this will be like coverage mesh,
      !    it will cover ocean target per process)
      !  map between rof 2 ocn is in  mbrmapro ;
      ! after this, the sending of tags for second hop (ocn context) will use the new par comm graph,
      !  that has more precise info, that got created
         call seq_comm_getData(CPLID,  mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

         call seq_comm_getData(CPLID ,mpigrp=mpigrp_CPLID)   !  second group, the coupler group CPLID is global variable

         type1 = 3 ! fv mesh nowadays
         direction = 1 !
         context_id = ocn(1)%cplcompid
         ! this creates a par comm graph between mbrxid and mbrxoid, with ids rof(1)%cplcompid, context ocn(1)%cplcompid
         ! this will be used in send/receive mappers
         ierr = iMOAB_MigrateMapMesh (mbrxid, mbrmapro, mbrxoid, mpicom_CPLID, mpigrp_CPLID, &
            mpigrp_CPLID, type1, rof(1)%cplcompid, context_id, direction)

         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in migrating rof mesh for map rof c2 ocn '
            call shr_sys_abort(subname//' ERROR in migrating rof mesh for map rof c2 ocn ')
         endif
         if (iamroot_CPLID)  then
            write(logunit,*) subname,' migrated mesh for map rof 2 ocn '
         endif
         if (mbrxoid .ge. 0) then ! we are on coupler side pes
            tagname=trim(seq_flds_r2x_fields)//C_NULL_CHAR
            tagtype = 1 ! dense, double
            numco= 1 ! 1  scalar per node
            ierr = iMOAB_DefineTagStorage(mbrxoid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining ' // trim(seq_flds_r2x_fields) // ' tags on coupler side in MOAB'
               call shr_sys_abort(subname//' ERROR in defining MOAB tags ')
            endif
         endif

         if (mboxid .ge. 0) then ! we are on coupler side pes, for ocean mesh
            tagname=trim(seq_flds_r2x_fields)//C_NULL_CHAR
            tagtype = 1 ! dense, double
            numco= 1 ! 1  scalar per node
            ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining ' // trim(seq_flds_r2x_fields) // ' tags on coupler side in MOAB, for ocean app'
               call shr_sys_abort(subname//' ERROR in defining MOAB tags ')
            endif
         endif
         if (iamroot_CPLID)  then
            write(logunit,*) subname,' created moab tags for seq_flds_r2x_fields '
         endif

! find out the number of local elements in moab mesh ocean instance on coupler
         ierr  = iMOAB_GetMeshInfo ( mboxid, nvert, nvise, nbl, nsurf, nvisBC )
         if (ierr .ne. 0) then
            write(logunit,*) subname,' cant get size of ocn mesh'
            call shr_sys_abort(subname//' ERROR in getting size of ocn mesh')
         endif
         ! ocn is cell mesh on coupler side
         mlsize = nvise(1)
         ent_type = 1 ! cell
         ! zero out the values just for r2x fields, on ocean instance
         nrflds = mct_aVect_nRattr(r2x_ox(1)) ! this is the size of r2x_fields
         arrsize = nrflds*mlsize
         allocate (tmparray(arrsize)) ! mlsize is the size of local land
         ! do we need to zero out others or just river ?
         tmparray = 0._R8
         ierr = iMOAB_SetDoubleTagStorage(mboxid, tagname, arrsize , ent_type, tmparray)
         if (ierr .ne. 0) then
            write(logunit,*) subname,' cant zero out r2x tags on ocn'
            call shr_sys_abort(subname//' cant zero out r2x tags on ocn')
         endif
         deallocate (tmparray)


         ! now we have to populate the map with the right moab attributes, so that it does the right projection
#ifdef MOABDEBUG
         if (mbrxoid.ge.0) then  ! we are on coupler PEs
            call mpi_comm_rank(mpicom_CPLID, rank_on_cpl  , ierr)
            if (rank_on_cpl .lt. 4) then
               prefix_output = "rof_cov"//CHAR(0)
               ierr = iMOAB_WriteLocalMesh(mbrxoid, prefix_output)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in writing coverage mesh rof 2 ocn '
               endif
            endif
         endif
#endif
! now take care of the mapper for MOAB mapper_Rr2o_liq
            if ( mapper_Rr2o_liq%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Rr2o_liq%mbname) &
                             //' mapper_Rr2o_liq'
                endif
            endif
            mapper_Rr2o_liq%src_mbid = mbrxid
            mapper_Rr2o_liq%tgt_mbid = mbrxoid ! this is special, it will really need this coverage type mesh
            mapper_Rr2o_liq%intx_mbid = mbrmapro
            mapper_Rr2o_liq%src_context = rof(1)%cplcompid
            mapper_Rr2o_liq%intx_context = ocn(1)%cplcompid ! this context was used in migrate mesh
            wgtIdef = 'map-from-file'//C_NULL_CHAR
            mapper_Rr2o_liq%weight_identifier = wgtIdef
            mapper_Rr2o_liq%mbname = 'mapper_Rr2o_liq'
            mapper_Rr2o_liq%read_map = .true.
#endif

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rr2o_ice'
          end if
          ! is this the same map as above ?
          call seq_map_init_rcfile(mapper_Rr2o_ice, rof(1), ocn(1), &
               'seq_maps.rc', 'rof2ocn_ice_rmapname:', 'rof2ocn_ice_rmaptype:',samegrid_ro, &
               'mapper_Rr2o_ice  initialization',esmf_map_flag)
! us the same one for mapper_Rr2o_ice and mapper_Fr2o
#ifdef HAVE_MOAB
! now take care of the mapper for MOAB mapper_Rr2o_ice
            if ( mapper_Rr2o_ice%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Rr2o_ice%mbname) &
                             //' mapper_Rr2o_ice'
                endif
            endif
            mapper_Rr2o_ice%src_mbid = mbrxid
            mapper_Rr2o_ice%tgt_mbid = mbrxoid ! special
            mapper_Rr2o_ice%intx_mbid = mbrmapro
            mapper_Rr2o_ice%src_context = rof(1)%cplcompid
            mapper_Rr2o_ice%intx_context = ocn(1)%cplcompid ! this context was used in migrate mesh
            wgtIdef = 'map-from-file'//C_NULL_CHAR
            mapper_Rr2o_ice%weight_identifier = wgtIdef
            mapper_Rr2o_ice%mbname = 'mapper_Rr2o_ice'
            mapper_Rr2o_ice%read_map = .true.
#endif
          if (flood_present) then
             if (iamroot_CPLID) then
                write(logunit,*) ' '
                write(logunit,F00) 'Initializing mapper_Fr2o'
             end if
             call seq_map_init_rcfile( mapper_Fr2o, rof(1), ocn(1), &
                  'seq_maps.rc', 'rof2ocn_fmapname:', 'rof2ocn_fmaptype:',samegrid_ro, &
                  string='mapper_Fr2o initialization', esmf_map=esmf_map_flag)
#ifdef HAVE_MOAB
! now take care of the mapper for MOAB mapper_Fr2o
            if ( mapper_Fr2o%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_Fr2o%mbname) &
                             //' mapper_Fr2o'
                endif
            endif
               mapper_Fr2o%src_mbid = mbrxid
               mapper_Fr2o%tgt_mbid = mbrxoid ! special
               mapper_Fr2o%intx_mbid = mbrmapro
               mapper_Fr2o%src_context = rof(1)%cplcompid
               mapper_Fr2o%intx_context = ocn(1)%cplcompid ! this context was used in migrate mesh
               wgtIdef = 'map-from-file'//C_NULL_CHAR
               mapper_Fr2o%weight_identifier = wgtIdef
               mapper_Fr2o%mbname = 'mapper_Fr2o'
#endif
          endif
       endif
       call shr_sys_flush(logunit)

       if (glc_c2_ocn) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rg2o_liq'
          end if
          call seq_map_init_rcfile(mapper_Rg2o_liq, glc(1), ocn(1), &
               'seq_maps.rc', 'glc2ocn_liq_rmapname:', 'glc2ocn_liq_rmaptype:',samegrid_og, &
               'mapper_Rg2o_liq initialization',esmf_map_flag)

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rg2o_ice'
          end if
          call seq_map_init_rcfile(mapper_Rg2o_ice, glc(1), ocn(1), &
               'seq_maps.rc', 'glc2ocn_ice_rmapname:', 'glc2ocn_ice_rmaptype:',samegrid_og, &
               'mapper_Rg2o_ice initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)

       if (glcshelf_c2_ocn) then !ice shelf coupled properties
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sg2o'
          end if
          call seq_map_init_rcfile(mapper_Sg2o, glc(1), ocn(1), &
               'seq_maps.rc', 'glc2ocn_smapname:', 'glc2ocn_smaptype:',samegrid_og, &
               'mapper_Sg2o initialization',esmf_map_flag)
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fg2o'
          end if
          call seq_map_init_rcfile(mapper_Fg2o, glc(1), ocn(1), &
               'seq_maps.rc', 'glc2ocn_fmapname:', 'glc2ocn_fmaptype:',samegrid_og, &
               'mapper_Fg2o initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)

       if (wav_c2_ocn) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sw2o'
          end if
          call seq_map_init_rcfile(mapper_Sw2o, wav(1), ocn(1), &
               'seq_maps.rc', 'wav2ocn_smapname:', 'wav2ocn_smaptype:',samegrid_ow, &
               'mapper_Sw2o initialization')
       endif
       call shr_sys_flush(logunit)

    end if
  end subroutine prep_ocn_init

  !================================================================================================

  subroutine prep_ocn_accum(timer)
    !---------------------------------------------------------------
    ! Description
    ! Accumulate ocn inputs
    ! Form partial sum of tavg ocn inputs (virtual "send" to ocn)
    ! NOTE: this is done AFTER the call to the merge in prep_ocn_mrg
    !
    ! Arguments
    character(len=*)        , intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi
    type(mct_avect) , pointer   :: x2o_ox
    character(*)    , parameter :: subname = '(prep_ocn_accum)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer), barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       x2o_ox => component_get_x2c_cx(ocn(eoi))

       if (x2oacc_ox_cnt == 0) then
          call mct_avect_copy(x2o_ox, x2oacc_ox(eoi))
       else
          call mct_avect_accum(x2o_ox, x2oacc_ox(eoi))
       endif
    enddo
    x2oacc_ox_cnt = x2oacc_ox_cnt + 1
    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_accum

  subroutine prep_ocn_accum_moab()
    !---------------------------------------------------------------
    ! Description
    ! Accumulate ocn inputs
    ! Form partial sum of tavg ocn inputs (virtual "send" to ocn)
    ! NOTE: this is done AFTER the call to the merge in prep_ocn_mrg
    use iMOAB, only : iMOAB_GetDoubleTagStorage
    ! Arguments
    !
    ! Local Variables
    integer   :: ent_type, ierr
    character(CXX)  :: tagname
    character(*)    , parameter :: subname = '(prep_ocn_accum_moab)'
    !---------------------------------------------------------------

    ! this method is called after merge, so it is not really necessary, because
    ! x2o_om should be saved between these calls
    tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
    ent_type = 1  ! cell type
    ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrSize_x2o_om , ent_type, x2o_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting x2o_om array  ')
    endif


    if (x2oacc_om_cnt == 0) then
       x2oacc_om = x2o_om
    else
       x2oacc_om = x2oacc_om + x2o_om
    endif

    x2oacc_om_cnt = x2oacc_om_cnt + 1

  end subroutine prep_ocn_accum_moab


  !================================================================================================

  subroutine prep_ocn_accum_avg(timer_accum)
    !---------------------------------------------------------------
    ! Description
    ! Finish accumulation ocn inputs
    !
    ! Arguments
    character(len=*), intent(in)    :: timer_accum
    !
    ! Local Variables
    integer :: eoi
    type(mct_avect), pointer :: x2o_ox
    character(*), parameter  :: subname = '(prep_ocn_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_accum), barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       ! temporary formation of average
       if (x2oacc_ox_cnt > 1) then
          call mct_avect_avg(x2oacc_ox(eoi), x2oacc_ox_cnt)
       end if

       ! ***NOTE***THE FOLLOWING ACTUALLY MODIFIES x2o_ox
       x2o_ox   => component_get_x2c_cx(ocn(eoi))
       call mct_avect_copy(x2oacc_ox(eoi), x2o_ox)
    enddo
    x2oacc_ox_cnt = 0
    call t_drvstopf (trim(timer_accum))

  end subroutine prep_ocn_accum_avg

subroutine prep_ocn_accum_avg_moab()
    !---------------------------------------------------------------
    ! Description
    ! Finish accumulation ocn inputs
    !
    ! Arguments
    use iMOAB, only : iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
    ! Local Variables
    integer   :: ent_type, ierr
    integer noflds, lsize ! used for restart case only?
    character(CXX)  :: tagname
    character(*), parameter  :: subname = '(prep_ocn_accum_avg_moab)'
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
    !---------------------------------------------------------------

       ! temporary formation of average
       if (x2oacc_om_cnt > 1) then
          !call mct_avect_avg(x2oacc_ox(eoi), x2oacc_ox_cnt)
          x2oacc_om = 1./x2oacc_om_cnt * x2oacc_om
       end if

       if (.not. allocated(x2o_om)) then
          ! we could come here in the restart case; not sure why only for 
          ! the case ERS_Vmoab_T62_oQU120.CMPASO-NYF
          lsize = size(x2oacc_om, 1)
          noflds = size(x2oacc_om, 2)
          allocate (x2o_om(lsize, noflds))
          arrSize_x2o_om = noflds * lsize
          
       endif

       ! ***NOTE***THE FOLLOWING ACTUALLY MODIFIES x2o_om

       x2o_om   = x2oacc_om
       !call mct_avect_copy(x2oacc_ox(eoi), x2o_ox)
       ! modify the tags
       tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
       ent_type = 1  ! cell type
       ierr = iMOAB_SetDoubleTagStorage ( mboxid, tagname, arrSize_x2o_om , ent_type, x2o_om)
       if (ierr .ne. 0) then
            call shr_sys_abort(subname//' error in setting x2o_om array  ')
      endif
#ifdef MOABDEBUG
      if (mboxid .ge. 0 ) then !  we are on coupler pes, for sure
         write(lnum,"(I0.2)")num_moab_exports
         outfile = 'OcnCplAftAvg'//trim(lnum)//'.h5m'//C_NULL_CHAR
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR
         ierr = iMOAB_WriteMesh(mboxid, trim(outfile), trim(wopts))
      endif
#endif

    x2oacc_om_cnt = 0

  end subroutine prep_ocn_accum_avg_moab

  !================================================================================================

  subroutine prep_ocn_mrg(infodata, fractions_ox, xao_ox, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge all ocn inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_ox(:)
    type(mct_aVect)         , intent(in)    :: xao_ox(:) ! Atm-ocn fluxes, ocn grid, cpl pes
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eii, ewi, egi, eoi, eai, eri, exi, efi, emi
    real(R8)                 :: flux_epbalfact ! adjusted precip factor
    type(mct_avect), pointer :: x2o_ox
    integer                  :: cnt
    character(*), parameter  :: subname = '(prep_ocn_mrg)'
    !---------------------------------------------------------------

    call seq_infodata_GetData(infodata, &
         flux_epbalfact=flux_epbalfact)

    call t_drvstartf (trim(timer_mrg), barrier=mpicom_CPLID)

    ! Use emi here for instance averaging capability, num_inst_max = num_inst_ocn normally
    ! if NOT x2o_average, just fill each instance of component_get_x2c_cx(ocn(eoi))
    ! if     x2o_average, then computer merge into x2o_ox_inst and accumulate that to
    !                     component_get_x2c_cx(ocn(1)) and then average it at the end

    if (x2o_average) then
       x2o_ox   => component_get_x2c_cx(ocn(1))
       call mct_aVect_zero(x2o_ox)
    endif

    cnt = 0
    do emi = 1,num_inst_max
       ! Use fortran mod to address ensembles in merge
       eoi = mod((emi-1),num_inst_ocn) + 1
       eai = mod((emi-1),num_inst_atm) + 1
       eii = mod((emi-1),num_inst_ice) + 1
       eri = mod((emi-1),num_inst_rof) + 1
       ewi = mod((emi-1),num_inst_wav) + 1
       egi = mod((emi-1),num_inst_glc) + 1
       exi = mod((emi-1),num_inst_xao) + 1
       efi = mod((emi-1),num_inst_frc) + 1

       if (x2o_average) then
          x2o_ox   => x2o_ox_inst
       else
          x2o_ox   => component_get_x2c_cx(ocn(eoi))
       endif

       call prep_ocn_merge( flux_epbalfact, a2x_ox(eai), i2x_ox(eii), r2x_ox(eri),  &
            w2x_ox(ewi), g2x_ox(egi), xao_ox(exi), fractions_ox(efi), x2o_ox )

       if (x2o_average) then
          x2o_ox   => component_get_x2c_cx(ocn(1))
          call mct_aVect_accum(x2o_ox_inst, x2o_ox)
          cnt = cnt + 1
       endif
    enddo

    if (x2o_average) then
       x2o_ox   => component_get_x2c_cx(ocn(1))
       call mct_avect_avg(x2o_ox,cnt)
    endif

    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_ocn_mrg

subroutine prep_ocn_mrg_moab(infodata, xao_ox)

    use iMOAB , only : iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, &
     iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
    use seq_comm_mct , only : mboxid, mbofxid ! ocean and atm-ocean flux instances
    !---------------------------------------------------------------
    ! Description
    ! Merge all ocn inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)        , pointer , intent(in)    :: xao_ox(:) ! Atm-ocn fluxes, ocn grid, cpl pes; used here just for indexing

    ! temporary, to compile
    ! type(mct_aVect)      :: fractions_o

    type(mct_avect) , pointer   :: a2x_o  ! used just for indexing
    type(mct_avect) , pointer   :: i2x_o
    type(mct_avect) , pointer   :: r2x_o
    type(mct_avect) , pointer   :: x2o_o
    type(mct_aVect) , pointer   :: xao_o
    !---------------------------------------------------------------


    real(R8)                 :: flux_epbalfact ! adjusted precip factor

    ! will build x2o_om , similar to x2o_ox
    ! no averages, just one ocn instance
    ! start copy from prep_ocn_merge
     ! Local variables
    integer  :: n,ka,ki,ko,kr,kw,kx,kir,kor,i,i1,o1
    integer  :: kof,kif
    integer  :: lsize, arrsize ! for double arrays
    integer , save :: noflds,naflds,niflds,nrflds,nxflds!  ,ngflds,nwflds, no glacier or wave model
    real(R8) :: ifrac,ifracr
    real(R8) :: afrac,afracr
    real(R8) :: frac_sum
    real(R8) :: avsdr, anidr, avsdf, anidf   ! albedos
    real(R8) :: fswabsv, fswabsi             ! sw
    character(CL),allocatable :: field_ocn(:)   ! string converted to char
    character(CL),allocatable :: field_atm(:)   ! string converted to char
    character(CL),allocatable :: field_ice(:)   ! string converted to char
    character(CL),allocatable :: field_rof(:)   ! string converted to char
    !character(CL),allocatable :: field_wav(:)   ! string converted to char
    character(CL),allocatable :: field_xao(:)   ! string converted to char
    !character(CL),allocatable :: field_glc(:)   ! string converted to char
    character(CL),allocatable :: itemc_ocn(:)   ! string converted to char
    character(CL),allocatable :: itemc_atm(:)   ! string converted to char
    character(CL),allocatable :: itemc_ice(:)   ! string converted to char
    character(CL),allocatable :: itemc_rof(:)   ! string converted to char
    !character(CL),allocatable :: itemc_wav(:)   ! string converted to char
    character(CL),allocatable :: itemc_xao(:)   ! string converted to char
    !character(CL),allocatable :: itemc_g2x(:)   ! string converted to char
    integer, save :: index_a2x_Faxa_swvdr
    integer, save :: index_a2x_Faxa_swvdf
    integer, save :: index_a2x_Faxa_swndr
    integer, save :: index_a2x_Faxa_swndf
    integer, save :: index_i2x_Fioi_swpen
    integer, save :: index_xao_So_avsdr
    integer, save :: index_xao_So_anidr
    integer, save :: index_xao_So_avsdf
    integer, save :: index_xao_So_anidf
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_r2x_Forr_rofl
    integer, save :: index_r2x_Forr_rofi
    integer, save :: index_r2x_Forr_rofl_16O
    integer, save :: index_r2x_Forr_rofi_16O
    integer, save :: index_r2x_Forr_rofl_18O
    integer, save :: index_r2x_Forr_rofi_18O
    integer, save :: index_r2x_Forr_rofl_HDO
    integer, save :: index_r2x_Forr_rofi_HDO
    integer, save :: index_r2x_Flrr_flood
    integer, save :: index_g2x_Fogg_rofl
    integer, save :: index_g2x_Fogg_rofi
    integer, save :: index_x2o_Foxx_swnet
    integer, save :: index_x2o_Faxa_snow
    integer, save :: index_x2o_Faxa_rain
    integer, save :: index_x2o_Faxa_prec
    integer, save :: index_x2o_Foxx_rofl
    integer, save :: index_x2o_Foxx_rofi
    integer, save :: index_x2o_Sf_afrac
    integer, save :: index_x2o_Sf_afracr
    integer, save :: index_x2o_Foxx_swnet_afracr
    integer, save :: index_x2o_Foxx_rofl_16O
    integer, save :: index_x2o_Foxx_rofi_16O
    integer, save :: index_x2o_Foxx_rofl_18O
    integer, save :: index_x2o_Foxx_rofi_18O
    integer, save :: index_x2o_Foxx_rofl_HDO
    integer, save :: index_x2o_Foxx_rofi_HDO
    integer, save :: index_a2x_Faxa_snowc_16O
    integer, save :: index_a2x_Faxa_snowl_16O
    integer, save :: index_a2x_Faxa_rainc_16O
    integer, save :: index_a2x_Faxa_rainl_16O
    integer, save :: index_x2o_Faxa_rain_16O
    integer, save :: index_x2o_Faxa_snow_16O
    integer, save :: index_x2o_Faxa_prec_16O
    integer, save :: index_a2x_Faxa_snowc_18O
    integer, save :: index_a2x_Faxa_snowl_18O
    integer, save :: index_a2x_Faxa_rainc_18O
    integer, save :: index_a2x_Faxa_rainl_18O
    integer, save :: index_x2o_Faxa_rain_18O
    integer, save :: index_x2o_Faxa_snow_18O
    integer, save :: index_x2o_Faxa_prec_18O
    integer, save :: index_a2x_Faxa_snowc_HDO
    integer, save :: index_a2x_Faxa_snowl_HDO
    integer, save :: index_a2x_Faxa_rainc_HDO
    integer, save :: index_a2x_Faxa_rainl_HDO
    integer, save :: index_x2o_Faxa_rain_HDO
    integer, save :: index_x2o_Faxa_snow_HDO
    integer, save :: index_x2o_Faxa_prec_HDO


    logical :: iamroot
    logical, save, pointer :: amerge(:),imerge(:),xmerge(:)
    integer, save, pointer :: aindx(:), iindx(:), xindx(:)
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: i2x_sharedindices
    type(mct_aVect_sharedindices),save :: r2x_sharedindices
    type(mct_aVect_sharedindices),save :: w2x_sharedindices
    type(mct_aVect_sharedindices),save :: xao_sharedindices
    type(mct_aVect_sharedindices),save :: g2x_sharedindices
    logical, save :: first_time = .true.

    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info

    character(CXX) ::tagname
    integer :: ent_type, ierr
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
#ifdef MOABCOMP
    character(CXX) :: mct_field
    real(R8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
#endif

! for moab, local allocatable arrays for each field, size of local ocean mesh
! these are the fields that are merged, in general
! some fields are already on the ocean instance (coming from projection)
!  (usually those on shared indices )
! all the rest will be needed for computation
! arrays will be allocated the first time, then filled with get tag values, merged, and set back to x2o ocean fields

    character(*),parameter :: subName = '(prep_ocn_merge_moab) '
    !-----------------------------------------------------------------------
    ! fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
    kif = 2 ! kif = mct_aVect_indexRa(fractions_o,"ifrac",perrWith=subName)
    kof = 3  !  kof = mct_aVect_indexRa(fractions_o,"ofrac",perrWith=subName)
    kir = 4  !  kir = mct_aVect_indexRa(fractions_o,"ifrad",perrWith=subName)
    kor = 5  !  kor = mct_aVect_indexRa(fractions_o,"ofrad",perrWith=subName)
    call seq_infodata_GetData(infodata, &
         flux_epbalfact=flux_epbalfact)

    call seq_comm_setptrs(CPLID, iamroot=iamroot)

 ! find out the number of local elements in moab mesh ocean instance on coupler
    ierr  = iMOAB_GetMeshInfo ( mboxid, nvert, nvise, nbl, nsurf, nvisBC )
    if (ierr .ne. 0) then
         write(logunit,*) subname,' error in getting info '
         call shr_sys_abort(subname//' error in getting info ')
    endif
    lsize = nvise(1) ! number of active cells

    if (first_time) then

       ! mct avs are used just for their fields metadata, not the actual reals
       ! (name of the fields)
       ! need these always, not only the first time
      a2x_o => a2x_ox(1)
      i2x_o => i2x_ox(1)
      r2x_o => r2x_ox(1)
      xao_o => xao_ox(1)
      x2o_o => component_get_x2c_cx(ocn(1))
      noflds = mct_aVect_nRattr(x2o_o) ! these are saved after first time
      naflds = mct_aVect_nRattr(a2x_o)
      niflds = mct_aVect_nRattr(i2x_o)
      nrflds = mct_aVect_nRattr(r2x_o)
   !nwflds = mct_aVect_nRattr(w2x_o)
      nxflds = mct_aVect_nRattr(xao_o)

      if (.not. allocated(x2o_om)) then
         !ngflds = mct_aVect_nRattr(g2x_o)
         allocate(x2o_om (lsize, noflds))
         arrSize_x2o_om = lsize * noflds ! this willbe used to set/get x2o_om tags
      endif
       allocate(a2x_om (lsize, naflds))
       allocate(i2x_om (lsize, niflds))
       allocate(r2x_om (lsize, nrflds))
       r2x_om = 0._R8 ! should we zero out all of them ?
       allocate(xao_om (lsize, nxflds))
       ! allocate fractions too
       ! use the fraclist fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
       allocate(fractions_om(lsize,5)) ! there are 5 fields here

       index_a2x_Faxa_swvdr     = mct_aVect_indexRA(a2x_o,'Faxa_swvdr')
       index_a2x_Faxa_swvdf     = mct_aVect_indexRA(a2x_o,'Faxa_swvdf')
       index_a2x_Faxa_swndr     = mct_aVect_indexRA(a2x_o,'Faxa_swndr')
       index_a2x_Faxa_swndf     = mct_aVect_indexRA(a2x_o,'Faxa_swndf')
       index_i2x_Fioi_swpen     = mct_aVect_indexRA(i2x_o,'Fioi_swpen')
       index_xao_So_avsdr       = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_anidr       = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_avsdf       = mct_aVect_indexRA(xao_o,'So_avsdf')
       index_xao_So_anidf       = mct_aVect_indexRA(xao_o,'So_anidf')
       index_x2o_Foxx_swnet     = mct_aVect_indexRA(x2o_o,'Foxx_swnet')

       index_a2x_Faxa_snowc     = mct_aVect_indexRA(a2x_o,'Faxa_snowc')
       index_a2x_Faxa_snowl     = mct_aVect_indexRA(a2x_o,'Faxa_snowl')
       index_a2x_Faxa_rainc     = mct_aVect_indexRA(a2x_o,'Faxa_rainc')
       index_a2x_Faxa_rainl     = mct_aVect_indexRA(a2x_o,'Faxa_rainl')
       index_r2x_Forr_rofl      = mct_aVect_indexRA(r2x_o,'Forr_rofl')
       index_r2x_Forr_rofi      = mct_aVect_indexRA(r2x_o,'Forr_rofi')
       index_r2x_Flrr_flood     = mct_aVect_indexRA(r2x_o,'Flrr_flood')
       !index_g2x_Fogg_rofl      = mct_aVect_indexRA(g2x_o,'Fogg_rofl')
       !index_g2x_Fogg_rofi      = mct_aVect_indexRA(g2x_o,'Fogg_rofi')
       index_x2o_Faxa_snow      = mct_aVect_indexRA(x2o_o,'Faxa_snow')
       index_x2o_Faxa_rain      = mct_aVect_indexRA(x2o_o,'Faxa_rain')
       index_x2o_Faxa_prec      = mct_aVect_indexRA(x2o_o,'Faxa_prec')
       index_x2o_Foxx_rofl      = mct_aVect_indexRA(x2o_o,'Foxx_rofl')
       index_x2o_Foxx_rofi      = mct_aVect_indexRA(x2o_o,'Foxx_rofi')

       if (seq_flds_i2o_per_cat) then
          index_x2o_Sf_afrac          = mct_aVect_indexRA(x2o_o,'Sf_afrac')
          index_x2o_Sf_afracr         = mct_aVect_indexRA(x2o_o,'Sf_afracr')
          index_x2o_Foxx_swnet_afracr = mct_aVect_indexRA(x2o_o,'Foxx_swnet_afracr')
       endif

       !wiso:
       ! H2_16O
       index_a2x_Faxa_snowc_16O = mct_aVect_indexRA(a2x_o,'Faxa_snowc_16O', perrWith='quiet')
       index_a2x_Faxa_snowl_16O = mct_aVect_indexRA(a2x_o,'Faxa_snowl_16O', perrWith='quiet')
       index_a2x_Faxa_rainc_16O = mct_aVect_indexRA(a2x_o,'Faxa_rainc_16O', perrWith='quiet')
       index_a2x_Faxa_rainl_16O = mct_aVect_indexRA(a2x_o,'Faxa_rainl_16O', perrWith='quiet')
       index_r2x_Forr_rofl_16O  = mct_aVect_indexRA(r2x_o,'Forr_rofl_16O' , perrWith='quiet')
       index_r2x_Forr_rofi_16O  = mct_aVect_indexRA(r2x_o,'Forr_rofi_16O' , perrWith='quiet')
       index_x2o_Faxa_rain_16O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_16O' , perrWith='quiet')
       index_x2o_Faxa_snow_16O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_16O' , perrWith='quiet')
       index_x2o_Faxa_prec_16O  = mct_aVect_indexRA(x2o_o,'Faxa_prec_16O' , perrWith='quiet')
       index_x2o_Foxx_rofl_16O  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_16O' , perrWith='quiet')
       index_x2o_Foxx_rofi_16O  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_16O' , perrWith='quiet')
       ! H2_18O
       index_a2x_Faxa_snowc_18O = mct_aVect_indexRA(a2x_o,'Faxa_snowc_18O', perrWith='quiet')
       index_a2x_Faxa_snowl_18O = mct_aVect_indexRA(a2x_o,'Faxa_snowl_18O', perrWith='quiet')
       index_a2x_Faxa_rainc_18O = mct_aVect_indexRA(a2x_o,'Faxa_rainc_18O', perrWith='quiet')
       index_a2x_Faxa_rainl_18O = mct_aVect_indexRA(a2x_o,'Faxa_rainl_18O', perrWith='quiet')
       index_r2x_Forr_rofl_18O  = mct_aVect_indexRA(r2x_o,'Forr_rofl_18O' , perrWith='quiet')
       index_r2x_Forr_rofi_18O  = mct_aVect_indexRA(r2x_o,'Forr_rofi_18O' , perrWith='quiet')
       index_x2o_Faxa_rain_18O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_18O' , perrWith='quiet')
       index_x2o_Faxa_snow_18O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_18O' , perrWith='quiet')
       index_x2o_Faxa_prec_18O  = mct_aVect_indexRA(x2o_o,'Faxa_prec_18O' , perrWith='quiet')
       index_x2o_Foxx_rofl_18O  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_18O' , perrWith='quiet')
       index_x2o_Foxx_rofi_18O  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_18O' , perrWith='quiet')
       ! HDO
       index_a2x_Faxa_snowc_HDO = mct_aVect_indexRA(a2x_o,'Faxa_snowc_HDO', perrWith='quiet')
       index_a2x_Faxa_snowl_HDO = mct_aVect_indexRA(a2x_o,'Faxa_snowl_HDO', perrWith='quiet')
       index_a2x_Faxa_rainc_HDO = mct_aVect_indexRA(a2x_o,'Faxa_rainc_HDO', perrWith='quiet')
       index_a2x_Faxa_rainl_HDO = mct_aVect_indexRA(a2x_o,'Faxa_rainl_HDO', perrWith='quiet')
       index_r2x_Forr_rofl_HDO  = mct_aVect_indexRA(r2x_o,'Forr_rofl_HDO' , perrWith='quiet')
       index_r2x_Forr_rofi_HDO  = mct_aVect_indexRA(r2x_o,'Forr_rofi_HDO' , perrWith='quiet')
       index_x2o_Faxa_rain_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_rain_HDO' , perrWith='quiet')
       index_x2o_Faxa_snow_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_snow_HDO' , perrWith='quiet')
       index_x2o_Faxa_prec_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_prec_HDO' , perrWith='quiet')
       index_x2o_Foxx_rofl_HDO  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_HDO' , perrWith='quiet')
       index_x2o_Foxx_rofi_HDO  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_HDO' , perrWith='quiet')

       ! Compute all other quantities based on standardized naming convention (see below)
       ! Only ocn field states that have the name-prefix Sx_ will be merged
       ! Only field names have the same name-suffix (after the "_") will be merged
       !    (e.g. Si_fldname, Sa_fldname => merged to => Sx_fldname)
       ! All fluxes will be scaled by the corresponding afrac or ifrac
       !   EXCEPT for
       !    -- Faxa_snnet, Faxa_snow, Faxa_rain, Faxa_prec (derived)
       ! All i2x_o fluxes that have the name-suffix "Faii" (atm/ice fluxes) will be ignored
       ! - only ice fluxes that are Fioi_... will be used in the ocean merges

       allocate(aindx(noflds), amerge(noflds))
       allocate(iindx(noflds), imerge(noflds))
       allocate(xindx(noflds), xmerge(noflds))
       allocate(field_atm(naflds), itemc_atm(naflds))
       allocate(field_ice(niflds), itemc_ice(niflds))
       allocate(field_ocn(noflds), itemc_ocn(noflds))
       allocate(field_rof(nrflds), itemc_rof(nrflds))
       !allocate(field_wav(nwflds), itemc_wav(nwflds))
       allocate(field_xao(nxflds), itemc_xao(nxflds))
       !allocate(field_glc(ngflds), itemc_g2x(ngflds))
       allocate(mrgstr(noflds))
       aindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       amerge(:) = .true.
       imerge(:) = .true.
       xmerge(:) = .true.

       do ko = 1,noflds
          field_ocn(ko) = mct_aVect_getRList2c(ko, x2o_o)
          itemc_ocn(ko) = trim(field_ocn(ko)(scan(field_ocn(ko),'_'):))
       enddo
       do ka = 1,naflds
          field_atm(ka) = mct_aVect_getRList2c(ka, a2x_o)
          itemc_atm(ka) = trim(field_atm(ka)(scan(field_atm(ka),'_'):))
       enddo
       do ki = 1,niflds
          field_ice(ki) = mct_aVect_getRList2c(ki, i2x_o)
          itemc_ice(ki) = trim(field_ice(ki)(scan(field_ice(ki),'_'):))
       enddo
       do kr = 1,nrflds
          field_rof(kr) = mct_aVect_getRList2c(kr, r2x_o)
          itemc_rof(kr) = trim(field_rof(kr)(scan(field_rof(kr),'_'):))
       enddo
      !  do kw = 1,nwflds
      !     field_wav(kw) = mct_aVect_getRList2c(kw, w2x_o)
      !     itemc_wav(kw) = trim(field_wav(kw)(scan(field_wav(kw),'_'):))
      !  enddo
       do kx = 1,nxflds
          field_xao(kx) = mct_aVect_getRList2c(kx, xao_o)
          itemc_xao(kx) = trim(field_xao(kx)(scan(field_xao(kx),'_'):))
       enddo
      !  do kx = 1,ngflds
      !     field_glc(kx) = mct_aVect_getRList2c(kx, g2x_o)
      !     itemc_g2x(kx) = trim(field_glc(kx)(scan(field_glc(kx),'_'):))
      !  enddo

       call mct_aVect_setSharedIndices(a2x_o, x2o_o, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(i2x_o, x2o_o, i2x_SharedIndices)
       call mct_aVect_setSharedIndices(r2x_o, x2o_o, r2x_SharedIndices)
       !call mct_aVect_setSharedIndices(w2x_o, x2o_o, w2x_SharedIndices)
       call mct_aVect_setSharedIndices(xao_o, x2o_o, xao_SharedIndices)
       !call mct_aVect_setSharedIndices(g2x_o, x2o_o, g2x_SharedIndices)

       do ko = 1,noflds
          !--- document merge ---
          mrgstr(ko) = subname//'x2o%'//trim(field_ocn(ko))//' ='
          if (field_ocn(ko)(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through
          end if
          if (field_ocn(ko)(1:1) == 'S' .and. field_ocn(ko)(2:2) /= 'x') then
             cycle ! ignore all ocn states that do not have a Sx_ prefix
          end if
          if (trim(field_ocn(ko)) == 'Foxx_swnet' .or. &
               trim(field_ocn(ko)) == 'Faxa_snow'  .or. &
               trim(field_ocn(ko)) == 'Faxa_rain'  .or. &
               trim(field_ocn(ko)) == 'Faxa_prec'  )then
             cycle ! ignore swnet, snow, rain, prec - treated explicitly above
          end if
          if (index(field_ocn(ko), 'Faxa_snow_' ) == 1 .or. &
               index(field_ocn(ko), 'Faxa_rain_' ) == 1 .or. &
               index(field_ocn(ko), 'Faxa_prec_' ) == 1 )then
             cycle ! ignore isotope snow, rain, prec - treated explicitly above
          end if
          !          if (trim(field_ocn(ko)(1:5)) == 'Foxx_') then
          !             cycle ! ignore runoff fields from land - treated in coupler
          !          end if

          do ka = 1,naflds
             if (trim(itemc_ocn(ko)) == trim(itemc_atm(ka))) then
                if ((trim(field_ocn(ko)) == trim(field_atm(ka)))) then
                   if (field_atm(ka)(1:1) == 'F') amerge(ko) = .false.
                end if
                ! --- make sure only one field matches ---
                if (aindx(ko) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ka field matches for ',trim(itemc_atm(ka))
                   call shr_sys_abort(subname//' ERROR multiple ka field matches')
                endif
                aindx(ko) = ka
             end if
          end do
          do ki = 1,niflds
             if (field_ice(ki)(1:1) == 'F' .and. field_ice(ki)(2:4) == 'aii') then
                cycle ! ignore all i2x_o fluxes that are ice/atm fluxes
             end if
             if (trim(itemc_ocn(ko)) == trim(itemc_ice(ki))) then
                if ((trim(field_ocn(ko)) == trim(field_ice(ki)))) then
                   if (field_ice(ki)(1:1) == 'F') imerge(ko) = .false.
                end if
                ! --- make sure only one field matches ---
                if (iindx(ko) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ki field matches for ',trim(itemc_ice(ki))
                   call shr_sys_abort(subname//' ERROR multiple ki field matches')
                endif
                iindx(ko) = ki
             end if
          end do
          do kx = 1,nxflds
             if (trim(itemc_ocn(ko)) == trim(itemc_xao(kx))) then
                if ((trim(field_ocn(ko)) == trim(field_xao(kx)))) then
                   if (field_xao(kx)(1:1) == 'F') xmerge(ko) = .false.
                end if
                ! --- make sure only one field matches ---
                if (xindx(ko) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple kx field matches for ',trim(itemc_xao(kx))
                   call shr_sys_abort(subname//' ERROR multiple kx field matches')
                endif
                xindx(ko) = kx
             end if
          end do

          ! --- add some checks ---

          ! --- make sure no merge of BOTH atm and xao ---
          if (aindx(ko) > 0 .and. xindx(ko) > 0) then
             write(logunit,*) subname,' ERROR: aindx and xindx both non-zero, not allowed'
             call shr_sys_abort(subname//' ERROR aindx and xindx both non-zero')
          endif

          ! --- make sure all terms agree on merge or non-merge aspect ---
          if (aindx(ko) > 0 .and. iindx(ko) > 0 .and. (amerge(ko) .neqv. imerge(ko))) then
             write(logunit,*) subname,' ERROR: aindx and iindx merge logic error'
             call shr_sys_abort(subname//' ERROR aindx and iindx merge logic error')
          endif
          if (aindx(ko) > 0 .and. xindx(ko) > 0 .and. (amerge(ko) .neqv. xmerge(ko))) then
             write(logunit,*) subname,' ERROR: aindx and xindx merge logic error'
             call shr_sys_abort(subname//' ERROR aindx and xindx merge logic error')
          endif
          if (xindx(ko) > 0 .and. iindx(ko) > 0 .and. (xmerge(ko) .neqv. imerge(ko))) then
             write(logunit,*) subname,' ERROR: xindx and iindx merge logic error'
             call shr_sys_abort(subname//' ERROR xindx and iindx merge logic error')
          endif

       end do

    end if

    !call mct_aVect_zero(x2o_o)
    ! replace with something else; make all x2o_fields 0 ? TODO

    !--- document copy operations ---
    if (first_time) then
       shared_fields_xao_x2o='' ! nothing in it yet
       !--- document merge ---
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field_atm(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A20, I3)" )i1, ' in a2x_o and x2o_o ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,i2x_SharedIndices%shared_real%num_indices
          i1=i2x_SharedIndices%shared_real%aVindices1(i)
          o1=i2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = i2x%'//trim(field_ice(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A20, I3)" )i1, ' in i2x_o and x2o_o ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,r2x_SharedIndices%shared_real%num_indices
          i1=r2x_SharedIndices%shared_real%aVindices1(i)
          o1=r2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = r2x%'//trim(field_rof(i1))
#ifdef MOABDEBUG
          write(lnum, "(I3, A20, I3)" )i1, ' in r2x_o and x2o_o ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
      !  do i=1,w2x_SharedIndices%shared_real%num_indices
      !     i1=w2x_SharedIndices%shared_real%aVindices1(i)
      !     o1=w2x_SharedIndices%shared_real%aVindices2(i)
      !     mrgstr(o1) = trim(mrgstr(o1))//' = w2x%'//trim(field_wav(i1))
      !  enddo
       do i=1,xao_SharedIndices%shared_real%num_indices
          i1=xao_SharedIndices%shared_real%aVindices1(i)
          o1=xao_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = xao%'//trim(field_xao(i1))
          ! will build tagname for moab set/get tag values
          shared_fields_xao_x2o = trim(shared_fields_xao_x2o)//trim(field_xao(i1))//':'
          size_of_shared_values = size_of_shared_values + lSize
#ifdef MOABDEBUG
          write(lnum, "(I3, A20, I3)" )i1, ' in xao_o and x2o_o ',o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       ! first time, allocate data for values_holder
       allocate(shared_values (size_of_shared_values))
      !  do i=1,g2x_SharedIndices%shared_real%num_indices
      !    i1=g2x_SharedIndices%shared_real%aVindices1(i)
      !    o1=g2x_SharedIndices%shared_real%aVindices2(i)
      !    mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field_glc(i1))
      ! enddo
    endif

    !    call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=r2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=w2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=mct_usevector)
    !call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    !call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=i2x_SharedIndices)
    !call mct_aVect_copy(aVin=r2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=r2x_SharedIndices)
    !!call mct_aVect_copy(aVin=w2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=w2x_SharedIndices)
    !call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=xao_SharedIndices)
    !!call mct_aVect_copy(aVin=g2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

    !--- document manual merges ---
    if (first_time) then
       mrgstr(index_x2o_Foxx_swnet) = trim(mrgstr(index_x2o_Foxx_swnet))//' = '// &
            'afracr*(a2x%Faxa_swvdr*(1.0-xao%So_avsdr) + '// &
            'a2x%Faxa_swvdf*(1.0-xao%So_avsdf) + '// &
            'a2x%Faxa_swndr*(1.0-xao%So_anidr) + '// &
            'a2x%Faxa_swndf*(1.0-xao%So_anidf)) + '// &
            'ifrac*i2x%Fioi_swpen'
       if (seq_flds_i2o_per_cat) then
          mrgstr(index_x2o_Foxx_swnet_afracr) = trim(mrgstr(index_x2o_Foxx_swnet_afracr))//' = '// &
               'afracr*(a2x%Faxa_swvdr*(1.0-xao%So_avsdr) + '// &
               'a2x%Faxa_swvdf*(1.0-xao%So_avsdf) + '// &
               'a2x%Faxa_swndr*(1.0-xao%So_anidr) + '// &
               'a2x%Faxa_swndf*(1.0-xao%So_anidf))'
       end if
       mrgstr(index_x2o_Faxa_snow) = trim(mrgstr(index_x2o_Faxa_snow))//' = '// &
            'afrac*(a2x%Faxa_snowc + a2x%Faxa_snowl)*flux_epbalfact'
       mrgstr(index_x2o_Faxa_rain) = trim(mrgstr(index_x2o_Faxa_rain))//' = '// &
            'afrac*(a2x%Faxa_rainc + a2x%Faxa_rainl)*flux_epbalfact'
       mrgstr(index_x2o_Faxa_prec) = trim(mrgstr(index_x2o_Faxa_prec))//' = '// &
            'afrac*(a2x%Faxa_snowc + a2x%Faxa_snowl + a2x%Faxa_rainc + a2x%Faxa_rainl)*flux_epbalfact'
       mrgstr(index_x2o_Foxx_rofl) = trim(mrgstr(index_x2o_Foxx_rofl))//' = '// &
            '(r2x%Forr_rofl + r2x%Flrr_flood + g2x%Fogg_rofl)*flux_epbalfact'
       mrgstr(index_x2o_Foxx_rofi) = trim(mrgstr(index_x2o_Foxx_rofi))//' = '// &
            '(r2x%Forr_rofi + g2x%Fogg_rofi)*flux_epbalfact'
       ! water isotope snow, rain prec
       if ( index_x2o_Faxa_snow_16O /= 0 )then
          mrgstr(index_x2o_Faxa_snow_16O) = trim(mrgstr(index_x2o_Faxa_snow_16O))//' = '// &
               'afrac*(a2x%Faxa_snowc_16O + a2x%Faxa_snowl_16O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_rain_16O) = trim(mrgstr(index_x2o_Faxa_rain_16O))//' = '// &
               'afrac*(a2x%Faxa_rainc_16O + a2x%Faxa_rainl_16O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_prec_16O) = trim(mrgstr(index_x2o_Faxa_prec_16O))//' = '// &
               'afrac*(a2x%Faxa_snowc_16O + a2x%Faxa_snowl_16O + a2x%Faxa_rainc_16O + '// &
               'a2x%Faxa_rainl_16O)*flux_epbalfact'
       end if
       if ( index_x2o_Faxa_snow_18O /= 0 )then
          mrgstr(index_x2o_Faxa_snow_18O) = trim(mrgstr(index_x2o_Faxa_snow_18O))//' = '// &
               'afrac*(a2x%Faxa_snowc_18O + a2x%Faxa_snowl_18O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_rain_18O) = trim(mrgstr(index_x2o_Faxa_rain_18O))//' = '// &
               'afrac*(a2x%Faxa_rainc_18O + a2x%Faxa_rainl_18O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_prec_18O) = trim(mrgstr(index_x2o_Faxa_prec_18O))//' = '// &
               'afrac*(a2x%Faxa_snowc_18O + a2x%Faxa_snowl_18O + a2x%Faxa_rainc_18O + '// &
               'a2x%Faxa_rainl_18O)*flux_epbalfact'
       end if
       if ( index_x2o_Faxa_snow_HDO /= 0 )then
          mrgstr(index_x2o_Faxa_snow_HDO) = trim(mrgstr(index_x2o_Faxa_snow_HDO))//' = '// &
               'afrac*(a2x%Faxa_snowc_HDO + a2x%Faxa_snowl_HDO)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_rain_HDO) = trim(mrgstr(index_x2o_Faxa_rain_HDO))//' = '// &
               'afrac*(a2x%Faxa_rainc_HDO + a2x%Faxa_rainl_HDO)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_prec_HDO) = trim(mrgstr(index_x2o_Faxa_prec_HDO))//' = '// &
               'afrac*(a2x%Faxa_snowc_HDO + a2x%Faxa_snowl_HDO + a2x%Faxa_rainc_HDO + '// &
               'a2x%Faxa_rainl_HDO)*flux_epbalfact'
       end if
    endif

    ! fill with fractions from ocean instance
    ent_type = 1 ! cells
    tagname = 'afrac:ifrac:ofrac:ifrad:ofrad'//C_NULL_CHAR
    arrsize = 5 * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, fractions_om)
    if (ierr .ne. 0) then
         call shr_sys_abort(subname//' error in getting fractions_om from ocean instance ')
    endif
   ! fill the o2x_om, etc double array fields noflds
    tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
    arrsize = noflds * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, x2o_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting x2o_om array ')
    endif
    ! zero out the output first (see line 1358)
    !x2o_om(:,:)=0.
    ! no, we should zero out only some indices, that accumulate
    do ko = 1, noflds
      if ( (aindx(ko) .gt. 0 ) .and. amerge(ko) ) then
         x2o_om(:, ko) = 0.
      endif
      if ( (iindx(ko) .gt. 0 ) .and. imerge(ko) ) then
         x2o_om(:, ko) = 0.
      endif
      if ( (xindx(ko) .gt. 0 ) .and. xmerge(ko) ) then
         x2o_om(:, ko) = 0.
      endif
    enddo
    tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
    arrsize = naflds * lsize !        allocate (a2x_om (lsize, naflds))
    ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, a2x_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting a2x_om array ')
    endif

    tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
    arrsize = niflds * lsize !        allocate (i2x_om (lsize, niflds))
    ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, i2x_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting i2x_om array ')
    endif

    tagname = trim(seq_flds_r2x_fields)//C_NULL_CHAR
    arrsize = nrflds * lsize !        allocate (r2x_om (lsize, nrflds))
    ierr = iMOAB_GetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, r2x_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting r2x_om array ')
    endif

    tagname = trim(seq_flds_xao_fields)//C_NULL_CHAR
    arrsize = nxflds * lsize !        allocate (xao_om (lsize, nxflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbofxid, tagname, arrsize , ent_type, xao_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting xao_om array ')
    endif

! #ifdef NOTDEF

    do n = 1,lsize

       ifrac = fractions_om(n,kif) ! fo_kif_ifrac(n) ! fractions_o%rAttr(kif,n)
       afrac = fractions_om(n,kof) ! fo_kof_ofrac(n) ! fractions_o%rAttr(kof,n)
       frac_sum = ifrac + afrac
       if ((frac_sum) /= 0._R8) then
          ifrac = ifrac / (frac_sum)
          afrac = afrac / (frac_sum)
       endif

       ifracr = fractions_om(n,kir) ! fo_kir_ifrad(n)  ! fractions_o%rAttr(kir,n)
       afracr = fractions_om(n,kor) ! fo_kor_ofrad(n) ! fractions_o%rAttr(kor,n)
       frac_sum = ifracr + afracr
       if ((frac_sum) /= 0._R8) then
          ifracr = ifracr / (frac_sum)
          afracr = afracr / (frac_sum)
       endif

       ! Derived: compute net short-wave
       avsdr = xao_om(n,index_xao_So_avsdr) ! avsdr = xao_So_avsdr(n) ! xao_o%rAttr(index_xao_So_avsdr,n)
       anidr = xao_om(n,index_xao_So_anidr) ! xao_So_anidr(n) !xao_o%rAttr(index_xao_So_anidr,n)
       avsdf = xao_om(n,index_xao_So_avsdf) !xao_So_avsdf(n) !xao_o%rAttr(index_xao_So_avsdf,n)
       anidf = xao_om(n,index_xao_So_anidf) ! xao_So_anidf(n)! xao_o%rAttr(index_xao_So_anidf,n)
       fswabsv  = a2x_om(n,index_a2x_Faxa_swvdr) * (1.0_R8 - avsdr) & ! a2x_Faxa_swvdr(n) * (1.0_R8 - avsdr) & !a2x_o%rAttr(index_a2x_Faxa_swvdr,n) * (1.0_R8 - avsdr) &
          + a2x_om(n,index_a2x_Faxa_swvdf) * (1.0_R8 - avsdf)!  + a2x_Faxa_swvdf(n) * (1.0_R8 - avsdf)  !+ a2x_o%rAttr(index_a2x_Faxa_swvdf,n) * (1.0_R8 - avsdf)
       fswabsi  = a2x_om(n,index_a2x_Faxa_swndr) * (1.0_R8 - anidr) & ! a2x_Faxa_swndr(n) * (1.0_R8 - anidr)  & !a2x_o%rAttr(index_a2x_Faxa_swndr,n) * (1.0_R8 - anidr) &
          + a2x_om(n,index_a2x_Faxa_swndf) * (1.0_R8 - anidf) !+ a2x_Faxa_swndf(n) * (1.0_R8 - anidf)  ! + a2x_o%rAttr(index_a2x_Faxa_swndf,n) * (1.0_R8 - anidf)
       x2o_om(n,index_x2o_Foxx_swnet) = (fswabsv + fswabsi) * afracr + & !x2o_Foxx_swnet(n) = (fswabsv + fswabsi)  * afracr + & !x2o_o%rAttr(index_x2o_Foxx_swnet,n) = (fswabsv + fswabsi)                 * afracr + &
          i2x_om(n,index_i2x_Fioi_swpen) * ifrac ! i2x_Fioi_swpen(n) * ifrac ! i2x_o%rAttr(index_i2x_Fioi_swpen,n) * ifrac

       if (seq_flds_i2o_per_cat) then
          x2o_om(n,index_x2o_Sf_afrac)          = afrac
          x2o_om(n,index_x2o_Sf_afracr)         = afracr
          x2o_om(n,index_x2o_Foxx_swnet_afracr) = (fswabsv + fswabsi)       * afracr
       end if

       ! Derived: compute total precipitation - scale total precip and runoff

       x2o_om(n,index_x2o_Faxa_snow ) = a2x_om(n,index_a2x_Faxa_snowc) * afrac + &
            a2x_om(n,index_a2x_Faxa_snowl) * afrac
       x2o_om(n,index_x2o_Faxa_rain ) = a2x_om(n,index_a2x_Faxa_rainc) * afrac + &
            a2x_om(n,index_a2x_Faxa_rainl) * afrac

       x2o_om(n,index_x2o_Faxa_snow ) = x2o_om(n,index_x2o_Faxa_snow ) * flux_epbalfact
       x2o_om(n,index_x2o_Faxa_rain ) = x2o_om(n,index_x2o_Faxa_rain ) * flux_epbalfact

       x2o_om(n,index_x2o_Faxa_prec ) = x2o_om(n,index_x2o_Faxa_rain ) + &
            x2o_om(n,index_x2o_Faxa_snow )

       x2o_om(n,index_x2o_Foxx_rofl) = (r2x_om(n,index_r2x_Forr_rofl ) + &
            r2x_om(n,index_r2x_Flrr_flood) )
           ! g2x_om(n,index_g2x_Fogg_rofl )) * flux_epbalfact
       x2o_om(n,index_x2o_Foxx_rofi) = (r2x_om(n,index_r2x_Forr_rofi ) ) * flux_epbalfact
          !  g2x_om(n,index_g2x_Fogg_rofi )) * flux_epbalfact


       if ( index_x2o_Foxx_rofl_16O /= 0 ) then
          x2o_om(n,index_x2o_Foxx_rofl_16O) = (r2x_om(n,index_r2x_Forr_rofl_16O) + &
               r2x_om(n,index_r2x_Flrr_flood) ) * flux_epbalfact
             !  g2x_om(n,index_g2x_Fogg_rofl )) * flux_epbalfact
          x2o_om(n,index_x2o_Foxx_rofi_16O) = (r2x_om(n,index_r2x_Forr_rofi_16O ) ) * flux_epbalfact
             !  g2x_om(n,index_g2x_Fogg_rofi )) * flux_epbalfact
          x2o_om(n,index_x2o_Foxx_rofl_18O) = (r2x_om(n,index_r2x_Forr_rofl_18O) + &
               r2x_om(n,index_r2x_Flrr_flood) ) * flux_epbalfact
              ! g2x_om(n,index_g2x_Fogg_rofl )) * flux_epbalfact
          x2o_om(n,index_x2o_Foxx_rofi_18O) = (r2x_om(n,index_r2x_Forr_rofi_18O ) ) * flux_epbalfact
               !g2x_om(n,index_g2x_Fogg_rofi )) * flux_epbalfact
          x2o_om(n,index_x2o_Foxx_rofl_HDO) = (r2x_om(n,index_r2x_Forr_rofl_HDO) + &
               r2x_om(n,index_r2x_Flrr_flood) ) * flux_epbalfact
               !g2x_om(n,index_g2x_Fogg_rofl )) * flux_epbalfact
          x2o_om(n,index_x2o_Foxx_rofi_HDO) = (r2x_om(n,index_r2x_Forr_rofi_HDO ) ) * flux_epbalfact
              ! g2x_om(n,index_g2x_Fogg_rofi )) * flux_epbalfact
       end if

       ! Derived: water isotopes total preciptiation and scaling

       if ( index_x2o_Faxa_snow_16O /= 0 )then
          x2o_om(n,index_x2o_Faxa_snow_16O ) = a2x_om(n,index_a2x_Faxa_snowc_16O) * afrac + &
               a2x_om(n,index_a2x_Faxa_snowl_16O) * afrac
          x2o_om(n,index_x2o_Faxa_rain_16O ) = a2x_om(n,index_a2x_Faxa_rainc_16O) * afrac + &
               a2x_om(n,index_a2x_Faxa_rainl_16O) * afrac

          x2o_om(n,index_x2o_Faxa_snow_16O ) = x2o_om(n,index_x2o_Faxa_snow_16O ) * flux_epbalfact
          x2o_om(n,index_x2o_Faxa_rain_16O ) = x2o_om(n,index_x2o_Faxa_rain_16O ) * flux_epbalfact

          x2o_om(n,index_x2o_Faxa_prec_16O ) = x2o_om(n,index_x2o_Faxa_rain_16O ) + &
               x2o_om(n,index_x2o_Faxa_snow_16O )
       end if

       if ( index_x2o_Faxa_snow_18O /= 0 )then
          x2o_om(n,index_x2o_Faxa_snow_18O ) = a2x_om(n,index_a2x_Faxa_snowc_18O) * afrac + &
               a2x_om(n,index_a2x_Faxa_snowl_18O) * afrac
          x2o_om(n,index_x2o_Faxa_rain_18O ) = a2x_om(n,index_a2x_Faxa_rainc_18O) * afrac + &
               a2x_om(n,index_a2x_Faxa_rainl_18O) * afrac

          x2o_om(n,index_x2o_Faxa_snow_18O ) = x2o_om(n,index_x2o_Faxa_snow_18O ) * flux_epbalfact
          x2o_om(n,index_x2o_Faxa_rain_18O ) = x2o_om(n,index_x2o_Faxa_rain_18O ) * flux_epbalfact

          x2o_om(n,index_x2o_Faxa_prec_18O ) = x2o_om(n,index_x2o_Faxa_rain_18O ) + &
               x2o_om(n,index_x2o_Faxa_snow_18O )
       end if

       if ( index_x2o_Faxa_snow_HDO /= 0 )then
          x2o_om(n,index_x2o_Faxa_snow_HDO ) = a2x_om(n,index_a2x_Faxa_snowc_HDO) * afrac + &
               a2x_om(n,index_a2x_Faxa_snowl_HDO) * afrac
          x2o_om(n,index_x2o_Faxa_rain_HDO ) = a2x_om(n,index_a2x_Faxa_rainc_HDO) * afrac + &
               a2x_om(n,index_a2x_Faxa_rainl_HDO) * afrac

          x2o_om(n,index_x2o_Faxa_snow_HDO ) = x2o_om(n,index_x2o_Faxa_snow_HDO ) * flux_epbalfact
          x2o_om(n,index_x2o_Faxa_rain_HDO ) = x2o_om(n,index_x2o_Faxa_rain_HDO ) * flux_epbalfact

          x2o_om(n,index_x2o_Faxa_prec_HDO ) = x2o_om(n,index_x2o_Faxa_rain_HDO ) + &
               x2o_om(n,index_x2o_Faxa_snow_HDO )
       end if
    end do
! #endif

    do ko = 1,noflds
       !--- document merge ---
       if (first_time) then
          if (iindx(ko) > 0) then
             if (imerge(ko)) then
                mrgstr(ko) = trim(mrgstr(ko))//' + ifrac*i2x%'//trim(field_ice(iindx(ko)))
             else
                mrgstr(ko) = trim(mrgstr(ko))//' = ifrac*i2x%'//trim(field_ice(iindx(ko)))
             end if
          end if
          if (aindx(ko) > 0) then
             if (amerge(ko)) then
                mrgstr(ko) = trim(mrgstr(ko))//' + afrac*a2x%'//trim(field_atm(aindx(ko)))
             else
                mrgstr(ko) = trim(mrgstr(ko))//' = afrac*a2x%'//trim(field_atm(aindx(ko)))
             end if
          end if
          if (xindx(ko) > 0) then
             if (xmerge(ko)) then
                mrgstr(ko) = trim(mrgstr(ko))//' + afrac*xao%'//trim(field_xao(xindx(ko)))
             else
                mrgstr(ko) = trim(mrgstr(ko))//' = afrac*xao%'//trim(field_xao(xindx(ko)))
             end if
          end if
       endif

!

       do n = 1,lsize
          ifrac = fractions_om(n,kif) !fo_kif_ifrac(n) ! fractions_o%rAttr(kif)
          afrac = fractions_om(n,kof) ! fo_kof_ofrac(n) ! fractions_o%rAttr(kof,n)
          frac_sum = ifrac + afrac
          if ((frac_sum) /= 0._R8) then
             ifrac = ifrac / (frac_sum)
             afrac = afrac / (frac_sum)
          endif
          if (iindx(ko) > 0) then
             if (imerge(ko)) then
                x2o_om(n, ko) = x2o_om(n, ko) + i2x_om(n, iindx(ko)) * ifrac
             else
                x2o_om(n, ko) = i2x_om(n,iindx(ko)) * ifrac
             end if
          end if
          if (aindx(ko) > 0) then
             if (amerge(ko)) then
                x2o_om(n,ko) = x2o_om(n,ko) + a2x_om(n,aindx(ko)) * afrac
             else
                x2o_om(n,ko) = a2x_om(n, aindx(ko)) * afrac
             end if
          end if
          if (xindx(ko) > 0) then
             if (xmerge(ko)) then
                x2o_om(n,ko) = x2o_om(n,ko) + xao_om(n,xindx(ko)) * afrac
             else
                x2o_om(n,ko) = xao_om(n,xindx(ko)) * afrac
             end if
          end if
       end do

    end do
! after we are done, set x2o_om to the mboxid

    tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
    arrsize = noflds * lsize
    ierr = iMOAB_SetDoubleTagStorage ( mboxid, tagname, arrsize , ent_type, x2o_om)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in setting x2o_om array ')
    endif

    ! we still need to get/set the shared fields between xao and x2o:
    tagname = trim(shared_fields_xao_x2o)//C_NULL_CHAR
    ierr = iMOAB_GetDoubleTagStorage ( mbofxid, tagname, size_of_shared_values , ent_type, shared_values)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting shared_values array ')
    endif
    ierr = iMOAB_SetDoubleTagStorage ( mboxid, tagname, size_of_shared_values , ent_type, shared_values)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in setting shared_values array on ocean instance')
    endif



#ifdef MOABCOMP
    x2o_o => component_get_x2c_cx(ocn(1))
    ! loop over all fields in seq_flds_x2o_fields
    call mct_list_init(temp_list ,seq_flds_x2o_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! cell for ocean
    if (iamroot) print *, num_moab_exports, trim(seq_flds_x2o_fields)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call compare_mct_av_moab_tag(ocn(1), x2o_o, mct_field,  mboxid, tagname, ent_type, difference, first_time)
    enddo
    call mct_list_clean(temp_list)
#endif

#ifdef MOABDEBUG
    if (mboxid .ge. 0 ) then !  we are on coupler pes, for sure
     write(lnum,"(I0.2)")num_moab_exports
     outfile = 'OcnCplAftMm'//trim(lnum)//'.h5m'//C_NULL_CHAR
     wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
     ierr = iMOAB_WriteMesh(mboxid, trim(outfile), trim(wopts))
     if (ierr .ne. 0) then
       call shr_sys_abort(subname//' error in writing ocean after merging')
     endif
   endif
#endif
    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do ko = 1,noflds
             write(logunit,'(A)') trim(mrgstr(ko))
          enddo
          write(logunit,'(A)') subname//' shared fields between xao and x2o '//trim(shared_fields_xao_x2o)
#ifdef MOABDEBUG
          write(logunit, *) ' Ocean fields computed on coupler'
          do ko=1,noflds
             write(logunit, *) trim(field_ocn(ko)), aindx(ko), amerge(ko), iindx(ko), imerge(ko), xindx(ko), xmerge(ko)
          enddo
          write(logunit, *) ' Atm fields projected on coupler'
          do ka = 1,naflds
             write(logunit, *) trim(field_atm(ka))
          enddo
          write(logunit, *) ' Ice fields projected on coupler'
          do ki = 1,niflds
             write(logunit, *) trim(field_ice(ki))
          enddo
          write(logunit, *) ' Runoff fields projected on coupler'
          do kr = 1,nrflds
             write(logunit, *) trim(field_rof(kr))
          enddo
          write(logunit, *) ' xao flux fields '
          do kx = 1,nxflds
             write(logunit, *) trim(field_xao(kx))
          enddo
#endif

       endif
       deallocate(mrgstr)

       deallocate(field_atm,itemc_atm)
       deallocate(field_ocn,itemc_ocn)
       deallocate(field_ice,itemc_ice)
       deallocate(field_rof,itemc_rof)
       !Sdeallocate(field_wav,itemc_wav)
       deallocate(field_xao,itemc_xao)
    endif

    first_time = .false.

    !end copy

  end subroutine prep_ocn_mrg_moab
  !================================================================================================

  subroutine prep_ocn_merge( flux_epbalfact, a2x_o, i2x_o, r2x_o, w2x_o, g2x_o, xao_o, &
       fractions_o, x2o_o )

    use prep_glc_mod, only: prep_glc_calculate_subshelf_boundary_fluxes

    !-----------------------------------------------------------------------
    !
    ! Arguments
    real(R8)       , intent(in)    :: flux_epbalfact
    type(mct_aVect), intent(in)    :: a2x_o
    type(mct_aVect), intent(in)    :: i2x_o
    type(mct_aVect), intent(in)    :: r2x_o
    type(mct_aVect), intent(in)    :: w2x_o
    type(mct_aVect), intent(in)    :: g2x_o
    type(mct_aVect), intent(in)    :: xao_o
    type(mct_aVect), intent(in)    :: fractions_o
    type(mct_aVect), intent(inout) :: x2o_o
    !
    ! Local variables
    integer  :: n,ka,ki,ko,kr,kw,kx,kir,kor,i,i1,o1
    integer  :: kof,kif
    integer  :: lsize
    integer  :: noflds,naflds,niflds,nrflds,nwflds,nxflds,ngflds
    real(R8) :: ifrac,ifracr
    real(R8) :: afrac,afracr
    real(R8) :: frac_sum
    real(R8) :: avsdr, anidr, avsdf, anidf   ! albedos
    real(R8) :: fswabsv, fswabsi             ! sw
    character(CL),allocatable :: field_ocn(:)   ! string converted to char
    character(CL),allocatable :: field_atm(:)   ! string converted to char
    character(CL),allocatable :: field_ice(:)   ! string converted to char
    character(CL),allocatable :: field_rof(:)   ! string converted to char
    character(CL),allocatable :: field_wav(:)   ! string converted to char
    character(CL),allocatable :: field_xao(:)   ! string converted to char
    character(CL),allocatable :: field_glc(:)   ! string converted to char
    character(CL),allocatable :: itemc_ocn(:)   ! string converted to char
    character(CL),allocatable :: itemc_atm(:)   ! string converted to char
    character(CL),allocatable :: itemc_ice(:)   ! string converted to char
    character(CL),allocatable :: itemc_rof(:)   ! string converted to char
    character(CL),allocatable :: itemc_wav(:)   ! string converted to char
    character(CL),allocatable :: itemc_xao(:)   ! string converted to char
    character(CL),allocatable :: itemc_g2x(:)   ! string converted to char
    integer, save :: index_a2x_Faxa_swvdr
    integer, save :: index_a2x_Faxa_swvdf
    integer, save :: index_a2x_Faxa_swndr
    integer, save :: index_a2x_Faxa_swndf
    integer, save :: index_i2x_Fioi_swpen
    integer, save :: index_xao_So_avsdr
    integer, save :: index_xao_So_anidr
    integer, save :: index_xao_So_avsdf
    integer, save :: index_xao_So_anidf
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_r2x_Forr_rofl
    integer, save :: index_r2x_Forr_rofi
    integer, save :: index_r2x_Forr_rofDIN
    integer, save :: index_r2x_Forr_rofDIP
    integer, save :: index_r2x_Forr_rofDON
    integer, save :: index_r2x_Forr_rofDOP
    integer, save :: index_r2x_Forr_rofDOC
    integer, save :: index_r2x_Forr_rofPP
    integer, save :: index_r2x_Forr_rofDSi
    integer, save :: index_r2x_Forr_rofPOC
    integer, save :: index_r2x_Forr_rofPN
    integer, save :: index_r2x_Forr_rofDIC
    integer, save :: index_r2x_Forr_rofFe
    integer, save :: index_r2x_Forr_rofl_16O
    integer, save :: index_r2x_Forr_rofi_16O
    integer, save :: index_r2x_Forr_rofl_18O
    integer, save :: index_r2x_Forr_rofi_18O
    integer, save :: index_r2x_Forr_rofl_HDO
    integer, save :: index_r2x_Forr_rofi_HDO
    integer, save :: index_r2x_Flrr_flood
    integer, save :: index_g2x_Fogg_rofl
    integer, save :: index_g2x_Fogg_rofi
    integer, save :: index_x2o_Foxx_swnet
    integer, save :: index_x2o_Faxa_snow
    integer, save :: index_x2o_Faxa_rain
    integer, save :: index_x2o_Faxa_prec
    integer, save :: index_x2o_Foxx_rofl
    integer, save :: index_x2o_Foxx_rofi
    integer, save :: index_x2o_Foxx_rofDIN
    integer, save :: index_x2o_Foxx_rofDIP
    integer, save :: index_x2o_Foxx_rofDON
    integer, save :: index_x2o_Foxx_rofDOP
    integer, save :: index_x2o_Foxx_rofDOC
    integer, save :: index_x2o_Foxx_rofPP
    integer, save :: index_x2o_Foxx_rofDSi
    integer, save :: index_x2o_Foxx_rofPOC
    integer, save :: index_x2o_Foxx_rofPN
    integer, save :: index_x2o_Foxx_rofDIC
    integer, save :: index_x2o_Foxx_rofFe
    integer, save :: index_x2o_Sf_afrac
    integer, save :: index_x2o_Sf_afracr
    integer, save :: index_x2o_Foxx_swnet_afracr
    integer, save :: index_x2o_Foxx_rofl_16O
    integer, save :: index_x2o_Foxx_rofi_16O
    integer, save :: index_x2o_Foxx_rofl_18O
    integer, save :: index_x2o_Foxx_rofi_18O
    integer, save :: index_x2o_Foxx_rofl_HDO
    integer, save :: index_x2o_Foxx_rofi_HDO
    integer, save :: index_a2x_Faxa_snowc_16O
    integer, save :: index_a2x_Faxa_snowl_16O
    integer, save :: index_a2x_Faxa_rainc_16O
    integer, save :: index_a2x_Faxa_rainl_16O
    integer, save :: index_x2o_Faxa_rain_16O
    integer, save :: index_x2o_Faxa_snow_16O
    integer, save :: index_x2o_Faxa_prec_16O
    integer, save :: index_a2x_Faxa_snowc_18O
    integer, save :: index_a2x_Faxa_snowl_18O
    integer, save :: index_a2x_Faxa_rainc_18O
    integer, save :: index_a2x_Faxa_rainl_18O
    integer, save :: index_x2o_Faxa_rain_18O
    integer, save :: index_x2o_Faxa_snow_18O
    integer, save :: index_x2o_Faxa_prec_18O
    integer, save :: index_a2x_Faxa_snowc_HDO
    integer, save :: index_a2x_Faxa_snowl_HDO
    integer, save :: index_a2x_Faxa_rainc_HDO
    integer, save :: index_a2x_Faxa_rainl_HDO
    integer, save :: index_x2o_Faxa_rain_HDO
    integer, save :: index_x2o_Faxa_snow_HDO
    integer, save :: index_x2o_Faxa_prec_HDO
    logical :: iamroot
    logical, save, pointer :: amerge(:),imerge(:),xmerge(:)
    integer, save, pointer :: aindx(:), iindx(:), xindx(:)
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: i2x_sharedindices
    type(mct_aVect_sharedindices),save :: r2x_sharedindices
    type(mct_aVect_sharedindices),save :: w2x_sharedindices
    type(mct_aVect_sharedindices),save :: xao_sharedindices
    type(mct_aVect_sharedindices),save :: g2x_sharedindices
    logical, save :: first_time = .true.
    character(*),parameter :: subName = '(prep_ocn_merge) '

    !-----------------------------------------------------------------------

    call seq_comm_setptrs(CPLID, iamroot=iamroot)

    noflds = mct_aVect_nRattr(x2o_o)
    naflds = mct_aVect_nRattr(a2x_o)
    niflds = mct_aVect_nRattr(i2x_o)
    nrflds = mct_aVect_nRattr(r2x_o)
    nwflds = mct_aVect_nRattr(w2x_o)
    nxflds = mct_aVect_nRattr(xao_o)
    ngflds = mct_aVect_nRattr(g2x_o)

    if (first_time) then
       index_a2x_Faxa_swvdr     = mct_aVect_indexRA(a2x_o,'Faxa_swvdr')
       index_a2x_Faxa_swvdf     = mct_aVect_indexRA(a2x_o,'Faxa_swvdf')
       index_a2x_Faxa_swndr     = mct_aVect_indexRA(a2x_o,'Faxa_swndr')
       index_a2x_Faxa_swndf     = mct_aVect_indexRA(a2x_o,'Faxa_swndf')
       index_i2x_Fioi_swpen     = mct_aVect_indexRA(i2x_o,'Fioi_swpen')
       index_xao_So_avsdr       = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_anidr       = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_avsdf       = mct_aVect_indexRA(xao_o,'So_avsdf')
       index_xao_So_anidf       = mct_aVect_indexRA(xao_o,'So_anidf')
       index_x2o_Foxx_swnet     = mct_aVect_indexRA(x2o_o,'Foxx_swnet')

       index_a2x_Faxa_snowc     = mct_aVect_indexRA(a2x_o,'Faxa_snowc')
       index_a2x_Faxa_snowl     = mct_aVect_indexRA(a2x_o,'Faxa_snowl')
       index_a2x_Faxa_rainc     = mct_aVect_indexRA(a2x_o,'Faxa_rainc')
       index_a2x_Faxa_rainl     = mct_aVect_indexRA(a2x_o,'Faxa_rainl')
       index_r2x_Forr_rofl      = mct_aVect_indexRA(r2x_o,'Forr_rofl')
       index_r2x_Forr_rofi      = mct_aVect_indexRA(r2x_o,'Forr_rofi')
       if (rof2ocn_nutrients) then
          index_r2x_Forr_rofDIN    = mct_aVect_indexRA(r2x_o,'Forr_rofDIN')
          index_r2x_Forr_rofDIP    = mct_aVect_indexRA(r2x_o,'Forr_rofDIP')
          index_r2x_Forr_rofDON    = mct_aVect_indexRA(r2x_o,'Forr_rofDON')
          index_r2x_Forr_rofDOP    = mct_aVect_indexRA(r2x_o,'Forr_rofDOP')
          index_r2x_Forr_rofDOC    = mct_aVect_indexRA(r2x_o,'Forr_rofDOC')
          index_r2x_Forr_rofPP     = mct_aVect_indexRA(r2x_o,'Forr_rofPP')
          index_r2x_Forr_rofDSi    = mct_aVect_indexRA(r2x_o,'Forr_rofDSi')
          index_r2x_Forr_rofPOC    = mct_aVect_indexRA(r2x_o,'Forr_rofPOC')
          index_r2x_Forr_rofPN     = mct_aVect_indexRA(r2x_o,'Forr_rofPN')
          index_r2x_Forr_rofDIC    = mct_aVect_indexRA(r2x_o,'Forr_rofDIC')
          index_r2x_Forr_rofFe     = mct_aVect_indexRA(r2x_o,'Forr_rofFe')
       endif
       index_r2x_Flrr_flood     = mct_aVect_indexRA(r2x_o,'Flrr_flood')
       index_g2x_Fogg_rofl      = mct_aVect_indexRA(g2x_o,'Fogg_rofl')
       index_g2x_Fogg_rofi      = mct_aVect_indexRA(g2x_o,'Fogg_rofi')
       index_x2o_Faxa_snow      = mct_aVect_indexRA(x2o_o,'Faxa_snow')
       index_x2o_Faxa_rain      = mct_aVect_indexRA(x2o_o,'Faxa_rain')
       index_x2o_Faxa_prec      = mct_aVect_indexRA(x2o_o,'Faxa_prec')
       index_x2o_Foxx_rofl      = mct_aVect_indexRA(x2o_o,'Foxx_rofl')
       index_x2o_Foxx_rofi      = mct_aVect_indexRA(x2o_o,'Foxx_rofi')
       if (rof2ocn_nutrients) then
          index_x2o_Foxx_rofDIN    = mct_aVect_indexRA(x2o_o,'Foxx_rofDIN')
          index_x2o_Foxx_rofDIP    = mct_aVect_indexRA(x2o_o,'Foxx_rofDIP')
          index_x2o_Foxx_rofDON    = mct_aVect_indexRA(x2o_o,'Foxx_rofDON')
          index_x2o_Foxx_rofDOP    = mct_aVect_indexRA(x2o_o,'Foxx_rofDOP')
          index_x2o_Foxx_rofDOC    = mct_aVect_indexRA(x2o_o,'Foxx_rofDOC')
          index_x2o_Foxx_rofPP     = mct_aVect_indexRA(x2o_o,'Foxx_rofPP')
          index_x2o_Foxx_rofDSi    = mct_aVect_indexRA(x2o_o,'Foxx_rofDSi')
          index_x2o_Foxx_rofPOC    = mct_aVect_indexRA(x2o_o,'Foxx_rofPOC')
          index_x2o_Foxx_rofPN     = mct_aVect_indexRA(x2o_o,'Foxx_rofPN')
          index_x2o_Foxx_rofDIC    = mct_aVect_indexRA(x2o_o,'Foxx_rofDIC')
          index_x2o_Foxx_rofFe     = mct_aVect_indexRA(x2o_o,'Foxx_rofFe')
       endif

       if (seq_flds_i2o_per_cat) then
          index_x2o_Sf_afrac          = mct_aVect_indexRA(x2o_o,'Sf_afrac')
          index_x2o_Sf_afracr         = mct_aVect_indexRA(x2o_o,'Sf_afracr')
          index_x2o_Foxx_swnet_afracr = mct_aVect_indexRA(x2o_o,'Foxx_swnet_afracr')
       endif

       !wiso:
       ! H2_16O
       index_a2x_Faxa_snowc_16O = mct_aVect_indexRA(a2x_o,'Faxa_snowc_16O', perrWith='quiet')
       index_a2x_Faxa_snowl_16O = mct_aVect_indexRA(a2x_o,'Faxa_snowl_16O', perrWith='quiet')
       index_a2x_Faxa_rainc_16O = mct_aVect_indexRA(a2x_o,'Faxa_rainc_16O', perrWith='quiet')
       index_a2x_Faxa_rainl_16O = mct_aVect_indexRA(a2x_o,'Faxa_rainl_16O', perrWith='quiet')
       index_r2x_Forr_rofl_16O  = mct_aVect_indexRA(r2x_o,'Forr_rofl_16O' , perrWith='quiet')
       index_r2x_Forr_rofi_16O  = mct_aVect_indexRA(r2x_o,'Forr_rofi_16O' , perrWith='quiet')
       index_x2o_Faxa_rain_16O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_16O' , perrWith='quiet')
       index_x2o_Faxa_snow_16O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_16O' , perrWith='quiet')
       index_x2o_Faxa_prec_16O  = mct_aVect_indexRA(x2o_o,'Faxa_prec_16O' , perrWith='quiet')
       index_x2o_Foxx_rofl_16O  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_16O' , perrWith='quiet')
       index_x2o_Foxx_rofi_16O  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_16O' , perrWith='quiet')
       ! H2_18O
       index_a2x_Faxa_snowc_18O = mct_aVect_indexRA(a2x_o,'Faxa_snowc_18O', perrWith='quiet')
       index_a2x_Faxa_snowl_18O = mct_aVect_indexRA(a2x_o,'Faxa_snowl_18O', perrWith='quiet')
       index_a2x_Faxa_rainc_18O = mct_aVect_indexRA(a2x_o,'Faxa_rainc_18O', perrWith='quiet')
       index_a2x_Faxa_rainl_18O = mct_aVect_indexRA(a2x_o,'Faxa_rainl_18O', perrWith='quiet')
       index_r2x_Forr_rofl_18O  = mct_aVect_indexRA(r2x_o,'Forr_rofl_18O' , perrWith='quiet')
       index_r2x_Forr_rofi_18O  = mct_aVect_indexRA(r2x_o,'Forr_rofi_18O' , perrWith='quiet')
       index_x2o_Faxa_rain_18O  = mct_aVect_indexRA(x2o_o,'Faxa_rain_18O' , perrWith='quiet')
       index_x2o_Faxa_snow_18O  = mct_aVect_indexRA(x2o_o,'Faxa_snow_18O' , perrWith='quiet')
       index_x2o_Faxa_prec_18O  = mct_aVect_indexRA(x2o_o,'Faxa_prec_18O' , perrWith='quiet')
       index_x2o_Foxx_rofl_18O  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_18O' , perrWith='quiet')
       index_x2o_Foxx_rofi_18O  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_18O' , perrWith='quiet')
       ! HDO
       index_a2x_Faxa_snowc_HDO = mct_aVect_indexRA(a2x_o,'Faxa_snowc_HDO', perrWith='quiet')
       index_a2x_Faxa_snowl_HDO = mct_aVect_indexRA(a2x_o,'Faxa_snowl_HDO', perrWith='quiet')
       index_a2x_Faxa_rainc_HDO = mct_aVect_indexRA(a2x_o,'Faxa_rainc_HDO', perrWith='quiet')
       index_a2x_Faxa_rainl_HDO = mct_aVect_indexRA(a2x_o,'Faxa_rainl_HDO', perrWith='quiet')
       index_r2x_Forr_rofl_HDO  = mct_aVect_indexRA(r2x_o,'Forr_rofl_HDO' , perrWith='quiet')
       index_r2x_Forr_rofi_HDO  = mct_aVect_indexRA(r2x_o,'Forr_rofi_HDO' , perrWith='quiet')
       index_x2o_Faxa_rain_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_rain_HDO' , perrWith='quiet')
       index_x2o_Faxa_snow_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_snow_HDO' , perrWith='quiet')
       index_x2o_Faxa_prec_HDO  = mct_aVect_indexRA(x2o_o,'Faxa_prec_HDO' , perrWith='quiet')
       index_x2o_Foxx_rofl_HDO  = mct_aVect_indexRA(x2o_o,'Foxx_rofl_HDO' , perrWith='quiet')
       index_x2o_Foxx_rofi_HDO  = mct_aVect_indexRA(x2o_o,'Foxx_rofi_HDO' , perrWith='quiet')

       ! Compute all other quantities based on standardized naming convention (see below)
       ! Only ocn field states that have the name-prefix Sx_ will be merged
       ! Only field names have the same name-suffix (after the "_") will be merged
       !    (e.g. Si_fldname, Sa_fldname => merged to => Sx_fldname)
       ! All fluxes will be scaled by the corresponding afrac or ifrac
       !   EXCEPT for
       !    -- Faxa_snnet, Faxa_snow, Faxa_rain, Faxa_prec (derived)
       ! All i2x_o fluxes that have the name-suffix "Faii" (atm/ice fluxes) will be ignored
       ! - only ice fluxes that are Fioi_... will be used in the ocean merges

       allocate(aindx(noflds), amerge(noflds))
       allocate(iindx(noflds), imerge(noflds))
       allocate(xindx(noflds), xmerge(noflds))
       allocate(field_atm(naflds), itemc_atm(naflds))
       allocate(field_ice(niflds), itemc_ice(niflds))
       allocate(field_ocn(noflds), itemc_ocn(noflds))
       allocate(field_rof(nrflds), itemc_rof(nrflds))
       allocate(field_wav(nwflds), itemc_wav(nwflds))
       allocate(field_xao(nxflds), itemc_xao(nxflds))
       allocate(field_glc(ngflds), itemc_g2x(ngflds))
       allocate(mrgstr(noflds))
       aindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       amerge(:) = .true.
       imerge(:) = .true.
       xmerge(:) = .true.

       do ko = 1,noflds
          field_ocn(ko) = mct_aVect_getRList2c(ko, x2o_o)
          itemc_ocn(ko) = trim(field_ocn(ko)(scan(field_ocn(ko),'_'):))
       enddo
       do ka = 1,naflds
          field_atm(ka) = mct_aVect_getRList2c(ka, a2x_o)
          itemc_atm(ka) = trim(field_atm(ka)(scan(field_atm(ka),'_'):))
       enddo
       do ki = 1,niflds
          field_ice(ki) = mct_aVect_getRList2c(ki, i2x_o)
          itemc_ice(ki) = trim(field_ice(ki)(scan(field_ice(ki),'_'):))
       enddo
       do kr = 1,nrflds
          field_rof(kr) = mct_aVect_getRList2c(kr, r2x_o)
          itemc_rof(kr) = trim(field_rof(kr)(scan(field_rof(kr),'_'):))
       enddo
       do kw = 1,nwflds
          field_wav(kw) = mct_aVect_getRList2c(kw, w2x_o)
          itemc_wav(kw) = trim(field_wav(kw)(scan(field_wav(kw),'_'):))
       enddo
       do kx = 1,nxflds
          field_xao(kx) = mct_aVect_getRList2c(kx, xao_o)
          itemc_xao(kx) = trim(field_xao(kx)(scan(field_xao(kx),'_'):))
       enddo
       do kx = 1,ngflds
          field_glc(kx) = mct_aVect_getRList2c(kx, g2x_o)
          itemc_g2x(kx) = trim(field_glc(kx)(scan(field_glc(kx),'_'):))
       enddo

       call mct_aVect_setSharedIndices(a2x_o, x2o_o, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(i2x_o, x2o_o, i2x_SharedIndices)
       call mct_aVect_setSharedIndices(r2x_o, x2o_o, r2x_SharedIndices)
       call mct_aVect_setSharedIndices(w2x_o, x2o_o, w2x_SharedIndices)
       call mct_aVect_setSharedIndices(xao_o, x2o_o, xao_SharedIndices)
       call mct_aVect_setSharedIndices(g2x_o, x2o_o, g2x_SharedIndices)

       do ko = 1,noflds
          !--- document merge ---
          mrgstr(ko) = subname//'x2o%'//trim(field_ocn(ko))//' ='
          if (field_ocn(ko)(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through
          end if
          if (field_ocn(ko)(1:1) == 'S' .and. field_ocn(ko)(2:2) /= 'x') then
             cycle ! ignore all ocn states that do not have a Sx_ prefix
          end if
          if (trim(field_ocn(ko)) == 'Foxx_swnet' .or. &
               trim(field_ocn(ko)) == 'Faxa_snow'  .or. &
               trim(field_ocn(ko)) == 'Faxa_rain'  .or. &
               trim(field_ocn(ko)) == 'Faxa_prec'  )then
             cycle ! ignore swnet, snow, rain, prec - treated explicitly above
          end if
          if (index(field_ocn(ko), 'Faxa_snow_' ) == 1 .or. &
               index(field_ocn(ko), 'Faxa_rain_' ) == 1 .or. &
               index(field_ocn(ko), 'Faxa_prec_' ) == 1 )then
             cycle ! ignore isotope snow, rain, prec - treated explicitly above
          end if
          !          if (trim(field_ocn(ko)(1:5)) == 'Foxx_') then
          !             cycle ! ignore runoff fields from land - treated in coupler
          !          end if

          do ka = 1,naflds
             if (trim(itemc_ocn(ko)) == trim(itemc_atm(ka))) then
                if ((trim(field_ocn(ko)) == trim(field_atm(ka)))) then
                   if (field_atm(ka)(1:1) == 'F') amerge(ko) = .false.
                end if
                ! --- make sure only one field matches ---
                if (aindx(ko) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ka field matches for ',trim(itemc_atm(ka))
                   call shr_sys_abort(subname//' ERROR multiple ka field matches')
                endif
                aindx(ko) = ka
             end if
          end do
          do ki = 1,niflds
             if (field_ice(ki)(1:1) == 'F' .and. field_ice(ki)(2:4) == 'aii') then
                cycle ! ignore all i2x_o fluxes that are ice/atm fluxes
             end if
             if (trim(itemc_ocn(ko)) == trim(itemc_ice(ki))) then
                if ((trim(field_ocn(ko)) == trim(field_ice(ki)))) then
                   if (field_ice(ki)(1:1) == 'F') imerge(ko) = .false.
                end if
                ! --- make sure only one field matches ---
                if (iindx(ko) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple ki field matches for ',trim(itemc_ice(ki))
                   call shr_sys_abort(subname//' ERROR multiple ki field matches')
                endif
                iindx(ko) = ki
             end if
          end do
          do kx = 1,nxflds
             if (trim(itemc_ocn(ko)) == trim(itemc_xao(kx))) then
                if ((trim(field_ocn(ko)) == trim(field_xao(kx)))) then
                   if (field_xao(kx)(1:1) == 'F') xmerge(ko) = .false.
                end if
                ! --- make sure only one field matches ---
                if (xindx(ko) /= 0) then
                   write(logunit,*) subname,' ERROR: found multiple kx field matches for ',trim(itemc_xao(kx))
                   call shr_sys_abort(subname//' ERROR multiple kx field matches')
                endif
                xindx(ko) = kx
             end if
          end do

          ! --- add some checks ---

          ! --- make sure no merge of BOTH atm and xao ---
          if (aindx(ko) > 0 .and. xindx(ko) > 0) then
             write(logunit,*) subname,' ERROR: aindx and xindx both non-zero, not allowed'
             call shr_sys_abort(subname//' ERROR aindx and xindx both non-zero')
          endif

          ! --- make sure all terms agree on merge or non-merge aspect ---
          if (aindx(ko) > 0 .and. iindx(ko) > 0 .and. (amerge(ko) .neqv. imerge(ko))) then
             write(logunit,*) subname,' ERROR: aindx and iindx merge logic error'
             call shr_sys_abort(subname//' ERROR aindx and iindx merge logic error')
          endif
          if (aindx(ko) > 0 .and. xindx(ko) > 0 .and. (amerge(ko) .neqv. xmerge(ko))) then
             write(logunit,*) subname,' ERROR: aindx and xindx merge logic error'
             call shr_sys_abort(subname//' ERROR aindx and xindx merge logic error')
          endif
          if (xindx(ko) > 0 .and. iindx(ko) > 0 .and. (xmerge(ko) .neqv. imerge(ko))) then
             write(logunit,*) subname,' ERROR: xindx and iindx merge logic error'
             call shr_sys_abort(subname//' ERROR xindx and iindx merge logic error')
          endif

       end do

    end if

    call mct_aVect_zero(x2o_o)

    !--- document copy operations ---
    if (first_time) then
       !--- document merge ---
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field_atm(i1))
       enddo
       do i=1,i2x_SharedIndices%shared_real%num_indices
          i1=i2x_SharedIndices%shared_real%aVindices1(i)
          o1=i2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = i2x%'//trim(field_ice(i1))
       enddo
       do i=1,r2x_SharedIndices%shared_real%num_indices
          i1=r2x_SharedIndices%shared_real%aVindices1(i)
          o1=r2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = r2x%'//trim(field_rof(i1))
       enddo
       do i=1,w2x_SharedIndices%shared_real%num_indices
          i1=w2x_SharedIndices%shared_real%aVindices1(i)
          o1=w2x_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = w2x%'//trim(field_wav(i1))
       enddo
       do i=1,xao_SharedIndices%shared_real%num_indices
          i1=xao_SharedIndices%shared_real%aVindices1(i)
          o1=xao_SharedIndices%shared_real%aVindices2(i)
          mrgstr(o1) = trim(mrgstr(o1))//' = xao%'//trim(field_xao(i1))
       enddo
       do i=1,g2x_SharedIndices%shared_real%num_indices
         i1=g2x_SharedIndices%shared_real%aVindices1(i)
         o1=g2x_SharedIndices%shared_real%aVindices2(i)
         mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field_glc(i1))
      enddo
    endif

    !    call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=r2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=w2x_o, aVout=x2o_o, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=i2x_SharedIndices)
    call mct_aVect_copy(aVin=r2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=r2x_SharedIndices)
    call mct_aVect_copy(aVin=w2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=w2x_SharedIndices)
    call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=xao_SharedIndices)
    call mct_aVect_copy(aVin=g2x_o, aVout=x2o_o, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

    !--- document manual merges ---
    if (first_time) then
       mrgstr(index_x2o_Foxx_swnet) = trim(mrgstr(index_x2o_Foxx_swnet))//' = '// &
            'afracr*(a2x%Faxa_swvdr*(1.0-xao%So_avsdr) + '// &
            'a2x%Faxa_swvdf*(1.0-xao%So_avsdf) + '// &
            'a2x%Faxa_swndr*(1.0-xao%So_anidr) + '// &
            'a2x%Faxa_swndf*(1.0-xao%So_anidf)) + '// &
            'ifrac*i2x%Fioi_swpen'
       if (seq_flds_i2o_per_cat) then
          mrgstr(index_x2o_Foxx_swnet_afracr) = trim(mrgstr(index_x2o_Foxx_swnet_afracr))//' = '// &
               'afracr*(a2x%Faxa_swvdr*(1.0-xao%So_avsdr) + '// &
               'a2x%Faxa_swvdf*(1.0-xao%So_avsdf) + '// &
               'a2x%Faxa_swndr*(1.0-xao%So_anidr) + '// &
               'a2x%Faxa_swndf*(1.0-xao%So_anidf))'
       end if
       mrgstr(index_x2o_Faxa_snow) = trim(mrgstr(index_x2o_Faxa_snow))//' = '// &
            'afrac*(a2x%Faxa_snowc + a2x%Faxa_snowl)*flux_epbalfact'
       mrgstr(index_x2o_Faxa_rain) = trim(mrgstr(index_x2o_Faxa_rain))//' = '// &
            'afrac*(a2x%Faxa_rainc + a2x%Faxa_rainl)*flux_epbalfact'
       mrgstr(index_x2o_Faxa_prec) = trim(mrgstr(index_x2o_Faxa_prec))//' = '// &
            'afrac*(a2x%Faxa_snowc + a2x%Faxa_snowl + a2x%Faxa_rainc + a2x%Faxa_rainl)*flux_epbalfact'
       mrgstr(index_x2o_Foxx_rofl) = trim(mrgstr(index_x2o_Foxx_rofl))//' = '// &
            '(r2x%Forr_rofl + r2x%Flrr_flood + g2x%Fogg_rofl)*flux_epbalfact'
       mrgstr(index_x2o_Foxx_rofi) = trim(mrgstr(index_x2o_Foxx_rofi))//' = '// &
            '(r2x%Forr_rofi + g2x%Fogg_rofi)*flux_epbalfact'

       ! river nutrients
       if (rof2ocn_nutrients) then
          mrgstr(index_x2o_Foxx_rofDIN) = trim(mrgstr(index_x2o_Foxx_rofDIN))//' = '// &
             'r2x%Forr_rofDIN'
          mrgstr(index_x2o_Foxx_rofDIP) = trim(mrgstr(index_x2o_Foxx_rofDIP))//' = '// &
             'r2x%Forr_rofDIP'
          mrgstr(index_x2o_Foxx_rofDON) = trim(mrgstr(index_x2o_Foxx_rofDON))//' = '// &
             'r2x%Forr_rofDON'
          mrgstr(index_x2o_Foxx_rofDOP) = trim(mrgstr(index_x2o_Foxx_rofDOP))//' = '// &
             'r2x%Forr_rofDOP'
          mrgstr(index_x2o_Foxx_rofDOC) = trim(mrgstr(index_x2o_Foxx_rofDOC))//' = '// &
             'r2x%Forr_rofDOC'
          mrgstr(index_x2o_Foxx_rofPP) = trim(mrgstr(index_x2o_Foxx_rofPP))//' = '// &
             'r2x%Forr_rofPP'
          mrgstr(index_x2o_Foxx_rofDSi) = trim(mrgstr(index_x2o_Foxx_rofDSi))//' = '// &
             'r2x%Forr_rofDSi'
          mrgstr(index_x2o_Foxx_rofPOC) = trim(mrgstr(index_x2o_Foxx_rofPOC))//' = '// &
             'r2x%Forr_rofPOC'
          mrgstr(index_x2o_Foxx_rofPN) = trim(mrgstr(index_x2o_Foxx_rofPN))//' = '// &
             'r2x%Forr_rofPN'
          mrgstr(index_x2o_Foxx_rofDIC) = trim(mrgstr(index_x2o_Foxx_rofDIC))//' = '// &
             'r2x%Forr_rofDIC'
          mrgstr(index_x2o_Foxx_rofFe) = trim(mrgstr(index_x2o_Foxx_rofFe))//' = '// &
             'r2x%Forr_rofFe'
       endif

       ! water isotope snow, rain prec
       if ( index_x2o_Faxa_snow_16O /= 0 )then
          mrgstr(index_x2o_Faxa_snow_16O) = trim(mrgstr(index_x2o_Faxa_snow_16O))//' = '// &
               'afrac*(a2x%Faxa_snowc_16O + a2x%Faxa_snowl_16O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_rain_16O) = trim(mrgstr(index_x2o_Faxa_rain_16O))//' = '// &
               'afrac*(a2x%Faxa_rainc_16O + a2x%Faxa_rainl_16O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_prec_16O) = trim(mrgstr(index_x2o_Faxa_prec_16O))//' = '// &
               'afrac*(a2x%Faxa_snowc_16O + a2x%Faxa_snowl_16O + a2x%Faxa_rainc_16O + '// &
               'a2x%Faxa_rainl_16O)*flux_epbalfact'
       end if
       if ( index_x2o_Faxa_snow_18O /= 0 )then
          mrgstr(index_x2o_Faxa_snow_18O) = trim(mrgstr(index_x2o_Faxa_snow_18O))//' = '// &
               'afrac*(a2x%Faxa_snowc_18O + a2x%Faxa_snowl_18O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_rain_18O) = trim(mrgstr(index_x2o_Faxa_rain_18O))//' = '// &
               'afrac*(a2x%Faxa_rainc_18O + a2x%Faxa_rainl_18O)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_prec_18O) = trim(mrgstr(index_x2o_Faxa_prec_18O))//' = '// &
               'afrac*(a2x%Faxa_snowc_18O + a2x%Faxa_snowl_18O + a2x%Faxa_rainc_18O + '// &
               'a2x%Faxa_rainl_18O)*flux_epbalfact'
       end if
       if ( index_x2o_Faxa_snow_HDO /= 0 )then
          mrgstr(index_x2o_Faxa_snow_HDO) = trim(mrgstr(index_x2o_Faxa_snow_HDO))//' = '// &
               'afrac*(a2x%Faxa_snowc_HDO + a2x%Faxa_snowl_HDO)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_rain_HDO) = trim(mrgstr(index_x2o_Faxa_rain_HDO))//' = '// &
               'afrac*(a2x%Faxa_rainc_HDO + a2x%Faxa_rainl_HDO)*flux_epbalfact'
          mrgstr(index_x2o_Faxa_prec_HDO) = trim(mrgstr(index_x2o_Faxa_prec_HDO))//' = '// &
               'afrac*(a2x%Faxa_snowc_HDO + a2x%Faxa_snowl_HDO + a2x%Faxa_rainc_HDO + '// &
               'a2x%Faxa_rainl_HDO)*flux_epbalfact'
       end if
    endif

    ! Compute input ocn state (note that this only applies to non-land portion of gridcell)

    kif = mct_aVect_indexRa(fractions_o,"ifrac",perrWith=subName)
    kof = mct_aVect_indexRa(fractions_o,"ofrac",perrWith=subName)
    kir = mct_aVect_indexRa(fractions_o,"ifrad",perrWith=subName)
    kor = mct_aVect_indexRa(fractions_o,"ofrad",perrWith=subName)
    lsize = mct_aVect_lsize(x2o_o)
    do n = 1,lsize

       ifrac = fractions_o%rAttr(kif,n)
       afrac = fractions_o%rAttr(kof,n)
       frac_sum = ifrac + afrac
       if ((frac_sum) /= 0._R8) then
          ifrac = ifrac / (frac_sum)
          afrac = afrac / (frac_sum)
       endif

       ifracr = fractions_o%rAttr(kir,n)
       afracr = fractions_o%rAttr(kor,n)
       frac_sum = ifracr + afracr
       if ((frac_sum) /= 0._R8) then
          ifracr = ifracr / (frac_sum)
          afracr = afracr / (frac_sum)
       endif

       ! Derived: compute net short-wave
       avsdr = xao_o%rAttr(index_xao_So_avsdr,n)
       anidr = xao_o%rAttr(index_xao_So_anidr,n)
       avsdf = xao_o%rAttr(index_xao_So_avsdf,n)
       anidf = xao_o%rAttr(index_xao_So_anidf,n)
       fswabsv  =  a2x_o%rAttr(index_a2x_Faxa_swvdr,n) * (1.0_R8 - avsdr) &
            + a2x_o%rAttr(index_a2x_Faxa_swvdf,n) * (1.0_R8 - avsdf)
       fswabsi  =  a2x_o%rAttr(index_a2x_Faxa_swndr,n) * (1.0_R8 - anidr) &
            + a2x_o%rAttr(index_a2x_Faxa_swndf,n) * (1.0_R8 - anidf)
       x2o_o%rAttr(index_x2o_Foxx_swnet,n) = (fswabsv + fswabsi)                 * afracr + &
            i2x_o%rAttr(index_i2x_Fioi_swpen,n) * ifrac

       if (seq_flds_i2o_per_cat) then
          x2o_o%rAttr(index_x2o_Sf_afrac,n)          = afrac
          x2o_o%rAttr(index_x2o_Sf_afracr,n)         = afracr
          x2o_o%rAttr(index_x2o_Foxx_swnet_afracr,n) = (fswabsv + fswabsi)       * afracr
       end if

       ! Derived: compute total precipitation - scale total precip and runoff

       x2o_o%rAttr(index_x2o_Faxa_snow ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc,n) * afrac + &
            a2x_o%rAttr(index_a2x_Faxa_snowl,n) * afrac
       x2o_o%rAttr(index_x2o_Faxa_rain ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc,n) * afrac + &
            a2x_o%rAttr(index_a2x_Faxa_rainl,n) * afrac

       x2o_o%rAttr(index_x2o_Faxa_snow ,n) = x2o_o%rAttr(index_x2o_Faxa_snow ,n) * flux_epbalfact
       x2o_o%rAttr(index_x2o_Faxa_rain ,n) = x2o_o%rAttr(index_x2o_Faxa_rain ,n) * flux_epbalfact

       x2o_o%rAttr(index_x2o_Faxa_prec ,n) = x2o_o%rAttr(index_x2o_Faxa_rain ,n) + &
            x2o_o%rAttr(index_x2o_Faxa_snow ,n)

       x2o_o%rAttr(index_x2o_Foxx_rofl, n) = (r2x_o%rAttr(index_r2x_Forr_rofl , n) + &
            r2x_o%rAttr(index_r2x_Flrr_flood, n) + &
            g2x_o%rAttr(index_g2x_Fogg_rofl , n)) * flux_epbalfact
       x2o_o%rAttr(index_x2o_Foxx_rofi, n) = (r2x_o%rAttr(index_r2x_Forr_rofi , n) + &
            g2x_o%rAttr(index_g2x_Fogg_rofi , n)) * flux_epbalfact

       if (rof2ocn_nutrients) then
          x2o_o%rAttr(index_x2o_Foxx_rofDIN, n) = r2x_o%rAttr(index_r2x_Forr_rofDIN , n)
          x2o_o%rAttr(index_x2o_Foxx_rofDIP, n) = r2x_o%rAttr(index_r2x_Forr_rofDIP , n)
          x2o_o%rAttr(index_x2o_Foxx_rofDON, n) = r2x_o%rAttr(index_r2x_Forr_rofDON , n)
          x2o_o%rAttr(index_x2o_Foxx_rofDOP, n) = r2x_o%rAttr(index_r2x_Forr_rofDOP , n)
          x2o_o%rAttr(index_x2o_Foxx_rofDOC, n) = r2x_o%rAttr(index_r2x_Forr_rofDOC , n)
          x2o_o%rAttr(index_x2o_Foxx_rofPP , n) = r2x_o%rAttr(index_r2x_Forr_rofPP  , n)
          x2o_o%rAttr(index_x2o_Foxx_rofDSi, n) = r2x_o%rAttr(index_r2x_Forr_rofDSi , n)
          x2o_o%rAttr(index_x2o_Foxx_rofPOC, n) = r2x_o%rAttr(index_r2x_Forr_rofPOC , n)
          x2o_o%rAttr(index_x2o_Foxx_rofPN , n) = r2x_o%rAttr(index_r2x_Forr_rofPN  , n)
          x2o_o%rAttr(index_x2o_Foxx_rofDIC, n) = r2x_o%rAttr(index_r2x_Forr_rofDIC , n)
          x2o_o%rAttr(index_x2o_Foxx_rofFe, n)  = r2x_o%rAttr(index_r2x_Forr_rofFe  , n)
       endif

       if ( index_x2o_Foxx_rofl_16O /= 0 ) then
          x2o_o%rAttr(index_x2o_Foxx_rofl_16O, n) = (r2x_o%rAttr(index_r2x_Forr_rofl_16O, n) + &
               r2x_o%rAttr(index_r2x_Flrr_flood, n) + &
               g2x_o%rAttr(index_g2x_Fogg_rofl , n)) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Foxx_rofi_16O, n) = (r2x_o%rAttr(index_r2x_Forr_rofi_16O , n) + &
               g2x_o%rAttr(index_g2x_Fogg_rofi , n)) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Foxx_rofl_18O, n) = (r2x_o%rAttr(index_r2x_Forr_rofl_18O, n) + &
               r2x_o%rAttr(index_r2x_Flrr_flood, n) + &
               g2x_o%rAttr(index_g2x_Fogg_rofl , n)) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Foxx_rofi_18O, n) = (r2x_o%rAttr(index_r2x_Forr_rofi_18O , n) + &
               g2x_o%rAttr(index_g2x_Fogg_rofi , n)) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Foxx_rofl_HDO, n) = (r2x_o%rAttr(index_r2x_Forr_rofl_HDO, n) + &
               r2x_o%rAttr(index_r2x_Flrr_flood, n) + &
               g2x_o%rAttr(index_g2x_Fogg_rofl , n)) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Foxx_rofi_HDO, n) = (r2x_o%rAttr(index_r2x_Forr_rofi_HDO , n) + &
               g2x_o%rAttr(index_g2x_Fogg_rofi , n)) * flux_epbalfact
       end if

       ! Derived: water isotopes total preciptiation and scaling

       if ( index_x2o_Faxa_snow_16O /= 0 )then
          x2o_o%rAttr(index_x2o_Faxa_snow_16O ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc_16O,n) * afrac + &
               a2x_o%rAttr(index_a2x_Faxa_snowl_16O,n) * afrac
          x2o_o%rAttr(index_x2o_Faxa_rain_16O ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc_16O,n) * afrac + &
               a2x_o%rAttr(index_a2x_Faxa_rainl_16O,n) * afrac

          x2o_o%rAttr(index_x2o_Faxa_snow_16O ,n) = x2o_o%rAttr(index_x2o_Faxa_snow_16O ,n) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Faxa_rain_16O ,n) = x2o_o%rAttr(index_x2o_Faxa_rain_16O ,n) * flux_epbalfact

          x2o_o%rAttr(index_x2o_Faxa_prec_16O ,n) = x2o_o%rAttr(index_x2o_Faxa_rain_16O ,n) + &
               x2o_o%rAttr(index_x2o_Faxa_snow_16O ,n)
       end if

       if ( index_x2o_Faxa_snow_18O /= 0 )then
          x2o_o%rAttr(index_x2o_Faxa_snow_18O ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc_18O,n) * afrac + &
               a2x_o%rAttr(index_a2x_Faxa_snowl_18O,n) * afrac
          x2o_o%rAttr(index_x2o_Faxa_rain_18O ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc_18O,n) * afrac + &
               a2x_o%rAttr(index_a2x_Faxa_rainl_18O,n) * afrac

          x2o_o%rAttr(index_x2o_Faxa_snow_18O ,n) = x2o_o%rAttr(index_x2o_Faxa_snow_18O ,n) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Faxa_rain_18O ,n) = x2o_o%rAttr(index_x2o_Faxa_rain_18O ,n) * flux_epbalfact

          x2o_o%rAttr(index_x2o_Faxa_prec_18O ,n) = x2o_o%rAttr(index_x2o_Faxa_rain_18O ,n) + &
               x2o_o%rAttr(index_x2o_Faxa_snow_18O ,n)
       end if

       if ( index_x2o_Faxa_snow_HDO /= 0 )then
          x2o_o%rAttr(index_x2o_Faxa_snow_HDO ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc_HDO,n) * afrac + &
               a2x_o%rAttr(index_a2x_Faxa_snowl_HDO,n) * afrac
          x2o_o%rAttr(index_x2o_Faxa_rain_HDO ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc_HDO,n) * afrac + &
               a2x_o%rAttr(index_a2x_Faxa_rainl_HDO,n) * afrac

          x2o_o%rAttr(index_x2o_Faxa_snow_HDO ,n) = x2o_o%rAttr(index_x2o_Faxa_snow_HDO ,n) * flux_epbalfact
          x2o_o%rAttr(index_x2o_Faxa_rain_HDO ,n) = x2o_o%rAttr(index_x2o_Faxa_rain_HDO ,n) * flux_epbalfact

          x2o_o%rAttr(index_x2o_Faxa_prec_HDO ,n) = x2o_o%rAttr(index_x2o_Faxa_rain_HDO ,n) + &
               x2o_o%rAttr(index_x2o_Faxa_snow_HDO ,n)
       end if
    end do

    do ko = 1,noflds
       !--- document merge ---
       if (first_time) then
          if (iindx(ko) > 0) then
             if (imerge(ko)) then
                mrgstr(ko) = trim(mrgstr(ko))//' + ifrac*i2x%'//trim(field_ice(iindx(ko)))
             else
                mrgstr(ko) = trim(mrgstr(ko))//' = ifrac*i2x%'//trim(field_ice(iindx(ko)))
             end if
          end if
          if (aindx(ko) > 0) then
             if (amerge(ko)) then
                mrgstr(ko) = trim(mrgstr(ko))//' + afrac*a2x%'//trim(field_atm(aindx(ko)))
             else
                mrgstr(ko) = trim(mrgstr(ko))//' = afrac*a2x%'//trim(field_atm(aindx(ko)))
             end if
          end if
          if (xindx(ko) > 0) then
             if (xmerge(ko)) then
                mrgstr(ko) = trim(mrgstr(ko))//' + afrac*xao%'//trim(field_xao(xindx(ko)))
             else
                mrgstr(ko) = trim(mrgstr(ko))//' = afrac*xao%'//trim(field_xao(xindx(ko)))
             end if
          end if
       endif

       do n = 1,lsize
          ifrac = fractions_o%rAttr(kif,n)
          afrac = fractions_o%rAttr(kof,n)
          frac_sum = ifrac + afrac
          if ((frac_sum) /= 0._R8) then
             ifrac = ifrac / (frac_sum)
             afrac = afrac / (frac_sum)
          endif
          if (iindx(ko) > 0) then
             if (imerge(ko)) then
                x2o_o%rAttr(ko,n) = x2o_o%rAttr(ko,n) + i2x_o%rAttr(iindx(ko),n) * ifrac
             else
                x2o_o%rAttr(ko,n) = i2x_o%rAttr(iindx(ko),n) * ifrac
             end if
          end if
          if (aindx(ko) > 0) then
             if (amerge(ko)) then
                x2o_o%rAttr(ko,n) = x2o_o%rAttr(ko,n) + a2x_o%rAttr(aindx(ko),n) * afrac
             else
                x2o_o%rAttr(ko,n) = a2x_o%rAttr(aindx(ko),n) * afrac
             end if
          end if
          if (xindx(ko) > 0) then
             if (xmerge(ko)) then
                x2o_o%rAttr(ko,n) = x2o_o%rAttr(ko,n) + xao_o%rAttr(xindx(ko),n) * afrac
             else
                x2o_o%rAttr(ko,n) = xao_o%rAttr(xindx(ko),n) * afrac
             end if
          end if
       end do
    end do

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do ko = 1,noflds
             write(logunit,'(A)') trim(mrgstr(ko))
          enddo
       endif
       deallocate(mrgstr)
       deallocate(field_atm,itemc_atm)
       deallocate(field_ocn,itemc_ocn)
       deallocate(field_ice,itemc_ice)
       deallocate(field_rof,itemc_rof)
       deallocate(field_wav,itemc_wav)
       deallocate(field_xao,itemc_xao)
    endif

    first_time = .false.

  end subroutine prep_ocn_merge

  !================================================================================================

  subroutine prep_ocn_calc_a2x_ox(timer)
    !---------------------------------------------------------------
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    type(mct_avect), pointer :: a2x_ax
    character(*), parameter  :: subname = '(prep_ocn_calc_a2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       a2x_ax => component_get_c2x_cx(atm(eai))

       call seq_map_map(mapper_Sa2o, a2x_ax, a2x_ox(eai), fldlist=seq_flds_a2x_states, norm=.true.)

       call seq_map_map(mapper_Fa2o, a2x_ax, a2x_ox(eai), fldlist=seq_flds_a2x_fluxes, norm=.true.)

#ifdef COMPARE_TO_NUOPC
       call seq_map_mapvect(mapper_Va2o, vect_map, a2x_ax, a2x_ox(eai), 'Sa_u', 'Sa_v', norm=.true.)
#else
       !--- tcx the norm should be true below, it's false for bfb backwards compatability
       call seq_map_mapvect(mapper_Va2o, vect_map, a2x_ax, a2x_ox(eai), 'Sa_u', 'Sa_v', norm=.false.)
#endif

    enddo

    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_calc_a2x_ox

  subroutine prep_ocn_calc_i2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create i2x_ox (note that i2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: eii
    type(mct_avect), pointer :: i2x_ix
    character(*), parameter  :: subname = '(prep_ocn_calc_i2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eii = 1,num_inst_ice
       i2x_ix => component_get_c2x_cx(ice(eii))
       call seq_map_map(mapper_SFi2o, i2x_ix, i2x_ox(eii), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_calc_i2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_r2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create r2x_ox (note that r2x_ox is a local module variable)
#ifdef MOABDEBUG   
    use iMOAB, only : iMOAB_WriteMesh
    use seq_comm_mct,        only: num_moab_exports  ! used to count the steps for moab files
#endif
    ! Arguments
    
    ! Local Variables
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
    integer         :: ierr
#endif
#ifdef MOABCOMP
     character*100             :: tagname, mct_field
     integer :: ent_type
     real*8 :: difference
#endif
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri
    type(mct_avect), pointer :: r2x_rx
    character(*), parameter  :: subname = '(prep_ocn_calc_r2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       r2x_rx => component_get_c2x_cx(rof(eri))
       call seq_map_map(mapper_Rr2o_liq, r2x_rx, r2x_ox(eri), &
            fldlist=seq_flds_r2o_liq_fluxes, norm=.false.)
       call seq_map_map(mapper_Rr2o_ice, r2x_rx, r2x_ox(eri), &
            fldlist=seq_flds_r2o_ice_fluxes, norm=.false.)
       if (flood_present) then
          call seq_map_map(mapper_Fr2o, r2x_rx, r2x_ox(eri), &
               fldlist='Flrr_flood', norm=.true.)
       endif
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_calc_r2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_g2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_ox (note that g2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_avect), pointer :: g2x_gx
    character(*),  parameter :: subname = '(prep_ocn_calc_g2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       g2x_gx => component_get_c2x_cx(glc(egi))
       call seq_map_map(mapper_Rg2o_liq, g2x_gx, g2x_ox(egi), &
            fldlist=seq_flds_g2o_liq_fluxes, norm=.false.)

       call seq_map_map(mapper_Rg2o_ice, g2x_gx, g2x_ox(egi), &
            fldlist=seq_flds_g2o_ice_fluxes, norm=.false.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_ocn_calc_g2x_ox

  !================================================================================================

  subroutine prep_ocn_shelf_calc_g2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_ox (note that g2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_avect), pointer :: g2x_gx
    character(*),  parameter :: subname = '(prep_ocn_calc_g2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       g2x_gx => component_get_c2x_cx(glc(egi))

       call seq_map_map(mapper_Sg2o, g2x_gx, g2x_ox(egi), norm=.true.)

       call seq_map_map(mapper_Fg2o, g2x_gx, g2x_ox(egi),norm=.true.)


    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_ocn_shelf_calc_g2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_w2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create w2x_ox (note that w2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: ewi
    type(mct_avect), pointer :: w2x_wx
    character(*), parameter  :: subname = '(prep_ocn_calc_w2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do ewi = 1,num_inst_wav
       w2x_wx => component_get_c2x_cx(wav(ewi))
       call seq_map_map(mapper_Sw2o, w2x_wx, w2x_ox(ewi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_ocn_calc_w2x_ox

  !================================================================================================

  function prep_ocn_get_a2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_a2x_ox(:)
    prep_ocn_get_a2x_ox => a2x_ox(:)
  end function prep_ocn_get_a2x_ox

  function prep_ocn_get_r2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_r2x_ox(:)
    prep_ocn_get_r2x_ox => r2x_ox(:)
  end function prep_ocn_get_r2x_ox

  function prep_ocn_get_i2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_i2x_ox(:)
    prep_ocn_get_i2x_ox => i2x_ox(:)
  end function prep_ocn_get_i2x_ox

  function prep_ocn_get_g2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_g2x_ox(:)
    prep_ocn_get_g2x_ox => g2x_ox(:)
  end function prep_ocn_get_g2x_ox

  function prep_ocn_get_w2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_w2x_ox(:)
    prep_ocn_get_w2x_ox => w2x_ox(:)
  end function prep_ocn_get_w2x_ox

  function prep_ocn_get_x2oacc_ox()
    type(mct_aVect), pointer :: prep_ocn_get_x2oacc_ox(:)
    prep_ocn_get_x2oacc_ox => x2oacc_ox(:)
  end function prep_ocn_get_x2oacc_ox

  function prep_ocn_get_x2oacc_ox_cnt()
    integer, pointer :: prep_ocn_get_x2oacc_ox_cnt
    prep_ocn_get_x2oacc_ox_cnt => x2oacc_ox_cnt
  end function prep_ocn_get_x2oacc_ox_cnt

  function prep_ocn_get_mapper_Sa2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Sa2o
    prep_ocn_get_mapper_Sa2o => mapper_Sa2o
  end function prep_ocn_get_mapper_Sa2o

  function prep_ocn_get_mapper_Va2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Va2o
    prep_ocn_get_mapper_Va2o => mapper_Va2o
  end function prep_ocn_get_mapper_Va2o

  function prep_ocn_get_mapper_Fa2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Fa2o
    prep_ocn_get_mapper_Fa2o => mapper_Fa2o
  end function prep_ocn_get_mapper_Fa2o

  function prep_ocn_get_mapper_Fr2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Fr2o
    prep_ocn_get_mapper_Fr2o => mapper_Fr2o
  end function prep_ocn_get_mapper_Fr2o

  function prep_ocn_get_mapper_Rr2o_liq()
    type(seq_map), pointer :: prep_ocn_get_mapper_Rr2o_liq
    prep_ocn_get_mapper_Rr2o_liq => mapper_Rr2o_liq
  end function prep_ocn_get_mapper_Rr2o_liq

  function prep_ocn_get_mapper_Rr2o_ice()
    type(seq_map), pointer :: prep_ocn_get_mapper_Rr2o_ice
    prep_ocn_get_mapper_Rr2o_ice => mapper_Rr2o_ice
  end function prep_ocn_get_mapper_Rr2o_ice

  function prep_ocn_get_mapper_SFi2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_SFi2o
    prep_ocn_get_mapper_SFi2o => mapper_SFi2o
  end function prep_ocn_get_mapper_SFi2o

  function prep_ocn_get_mapper_Rg2o_liq()
    type(seq_map), pointer :: prep_ocn_get_mapper_Rg2o_liq
    prep_ocn_get_mapper_Rg2o_liq => mapper_Rg2o_liq
  end function prep_ocn_get_mapper_Rg2o_liq

  function prep_ocn_get_mapper_Rg2o_ice()
    type(seq_map), pointer :: prep_ocn_get_mapper_Rg2o_ice
    prep_ocn_get_mapper_Rg2o_ice => mapper_Rg2o_ice
  end function prep_ocn_get_mapper_Rg2o_ice

  function prep_ocn_get_mapper_Sg2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Sg2o
    prep_ocn_get_mapper_Sg2o => mapper_Sg2o
  end function prep_ocn_get_mapper_Sg2o

  function prep_ocn_get_mapper_Fg2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Fg2o
    prep_ocn_get_mapper_Fg2o => mapper_Fg2o
  end function prep_ocn_get_mapper_Fg2o

  function prep_ocn_get_mapper_Sw2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Sw2o
    prep_ocn_get_mapper_Sw2o => mapper_Sw2o
  end function prep_ocn_get_mapper_Sw2o
  function prep_ocn_get_x2oacc_om()
    real(R8), DIMENSION(:, :), pointer :: prep_ocn_get_x2oacc_om
    prep_ocn_get_x2oacc_om => x2oacc_om
  end function prep_ocn_get_x2oacc_om
  function prep_ocn_get_x2oacc_om_cnt()
    integer, pointer :: prep_ocn_get_x2oacc_om_cnt
    prep_ocn_get_x2oacc_om_cnt => x2oacc_om_cnt
  end function prep_ocn_get_x2oacc_om_cnt

end module prep_ocn_mod
