module prep_ice_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_atm, num_inst_ocn, num_inst_glc
  use seq_comm_mct    , only: num_inst_ice, num_inst_frc, num_inst_rof
  use seq_comm_mct    , only: CPLID, ICEID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs
  use seq_comm_mct,     only: mboxid ! iMOAB id for mpas ocean migrated mesh to coupler pes
  use seq_comm_mct,     only : mbixid   ! iMOAB for sea-ice migrated to coupler
  use seq_comm_mct,     only : num_moab_exports

  use seq_comm_mct,     only : seq_comm_getinfo => seq_comm_setptrs

  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: ice, atm, ocn, glc, rof
  use iso_c_binding
#ifdef MOABCOMP
  use component_type_mod, only:  compare_mct_av_moab_tag
#endif

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_ice_init
  public :: prep_ice_mrg
  public :: prep_ice_mrg_moab

  public :: prep_ice_calc_a2x_ix
  public :: prep_ice_calc_o2x_ix
  public :: prep_ice_calc_r2x_ix
  public :: prep_ice_calc_g2x_ix
  public :: prep_ice_shelf_calc_g2x_ix

  public :: prep_ice_get_a2x_ix
  public :: prep_ice_get_o2x_ix
  public :: prep_ice_get_g2x_ix
  public :: prep_ice_get_r2x_ix

  public :: prep_ice_get_mapper_SFo2i
  public :: prep_ice_get_mapper_Rg2i
  public :: prep_ice_get_mapper_Sg2i
  public :: prep_ice_get_mapper_Fg2i

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_ice_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_SFo2i
  type(seq_map), pointer :: mapper_Rg2i
  type(seq_map), pointer :: mapper_Sg2i
  type(seq_map), pointer :: mapper_Fg2i
  type(seq_map), pointer :: mapper_Rr2i

  ! attribute vectors
  type(mct_aVect), pointer :: a2x_ix(:) ! Atm export, ice grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: o2x_ix(:) ! Ocn export, ice grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: g2x_ix(:) ! Glc export, ice grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: r2x_ix(:) ! Rof export, ice grid, cpl pes - allocated in driver

! MOAB arrays
  real (kind=r8) , allocatable, private :: x2i_im (:,:)
  real (kind=r8) , allocatable, private :: a2x_im (:,:)
  real (kind=r8) , allocatable, private :: r2x_im (:,:)
  ! seq_comm_getData variables
  integer :: mpicom_CPLID                         ! MPI cpl communicator
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_ice_init(infodata, ocn_c2_ice, glc_c2_ice, glcshelf_c2_ice, rof_c2_ice)

    use iMOAB, only: iMOAB_RegisterApplication, &
      iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_ComputeCommGraph
    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type (seq_infodata_type) , intent(in)    :: infodata
    logical,                   intent(in)    :: ocn_c2_ice ! .true.  => ocn to ice coupling on
    logical,                   intent(in)    :: glc_c2_ice ! .true.  => glc to ice coupling on
    logical,                   intent(in)    :: glcshelf_c2_ice ! .true.  => glc ice shelf to ice coupling on
    logical,                   intent(in)    :: rof_c2_ice ! .true.  => rof to ice coupling on
    !
    ! Local Variables
    integer                          :: lsize_i
    integer                          :: eai, eoi, egi, eri
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                          :: samegrid_ig   ! samegrid glc and ice
    logical                          :: samegrid_ro   ! samegrid rof and ice/ocn
    logical                          :: ice_present   ! .true. => ice is present
    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    character(CL)                    :: ice_gnam      ! ice grid
    character(CL)                    :: ocn_gnam      ! ocn grid
    character(CL)                    :: glc_gnam      ! glc grid
    character(CL)                    :: rof_gnam      ! rof grid
    type(mct_avect), pointer         :: i2x_ix
    character(*), parameter          :: subname = '(prep_ice_init)'
    logical                          :: no_match      ! to force a new map between ocean and ice, always
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
!MOAB stuff
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

    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         esmf_map_flag=esmf_map_flag  , &
         ice_present=ice_present      , &
         ice_gnam=ice_gnam            , &
         ocn_gnam=ocn_gnam            , &
         rof_gnam=rof_gnam            , &
         glc_gnam=glc_gnam)

    allocate(mapper_SFo2i)
    allocate(mapper_Rg2i)
    allocate(mapper_Sg2i)
    allocate(mapper_Fg2i)
    allocate(mapper_Rr2i)

    if (ice_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       i2x_ix => component_get_c2x_cx(ice(1))
       lsize_i = mct_aVect_lsize(i2x_ix)

       allocate(a2x_ix(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_ix(eai), rList=seq_flds_a2x_fields, lsize=lsize_i)
          call mct_aVect_zero(a2x_ix(eai))
       end do
       allocate(o2x_ix(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_ix(eoi), rList=seq_flds_o2x_fields, lsize=lsize_i)
          call mct_aVect_zero(o2x_ix(eoi))
       enddo
       allocate(g2x_ix(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(g2x_ix(egi), rList=seq_flds_g2x_fields, lsize=lsize_i)
          call mct_aVect_zero(g2x_ix(egi))
       enddo
       allocate(r2x_ix(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_aVect_init(r2x_ix(eri), rList=seq_flds_r2x_fields, lsize=lsize_i)
          call mct_aVect_zero(r2x_ix(eri))
       end do

       samegrid_ig = .true.
       samegrid_ro = .true.
       if (trim(ice_gnam) /= trim(glc_gnam)) samegrid_ig = .false.
       if (trim(rof_gnam) /= trim(ocn_gnam)) samegrid_ro = .false.

       if (ocn_c2_ice) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_SFo2i'
          end if
          no_match = .true.
          call seq_map_init_rearrolap(mapper_SFo2i, ocn(1), ice(1), 'mapper_SFo2i', no_match) ! force a new map always

#ifdef HAVE_MOAB
          if ( (mbixid .ge. 0) .and. (mboxid .ge. 0)) then
            ! moab also will do just a rearrange, hopefully, in this case, based on the comm graph
            !   that is computed here
            call seq_comm_getinfo(CPLID ,mpigrp=mpigrp_CPLID)   !  second group, the coupler group CPLID is global variable

            type1 = 3
            type2 = 3 ! fv-fv graph
            ! iMOAB compute comm graph ice-ocn, based on the same global id
            ! it will be a simple migrate from ice mesh directly to ocean, using the comm graph computed here
            ! TODO: find if CommGraph already exists.

            ierr = iMOAB_ComputeCommGraph( mboxid, mbixid, mpicom_CPLID, mpigrp_CPLID, mpigrp_CPLID, &
               type1, type2, ocn(1)%cplcompid, ice(1)%cplcompid)
            if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in computing graph ocn - ice x '
                  call shr_sys_abort(subname//' ERROR  in computing graph ocn -ice x ')
            endif

            ! define tags according to the seq_flds_i2x_fields
            tagtype = 1  ! dense, double
            numco = 1 !  one value per cell / entity
            tagname = trim(seq_flds_o2x_fields)//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
            if ( ierr == 1 ) then
               call shr_sys_abort( subname//' ERROR: cannot define tags for ice proj to ocn' )
            end if
            if ( mapper_SFo2i%src_mbid .gt. -1 ) then
                if (iamroot_CPLID) then
                     write(logunit,F00) 'overwriting '//trim(mapper_SFo2i%mbname) &
                             //' mapper_SFo2i'
                endif
            endif

            mapper_SFo2i%src_mbid = mboxid
            mapper_SFo2i%tgt_mbid = mbixid
            ! no intersection, so will have to transform data without it
            mapper_SFo2i%src_context = ocn(1)%cplcompid
            mapper_SFo2i%intx_context = ice(1)%cplcompid
            mapper_SFo2i%mbname = 'mapper_SFo2i'

            if(mapper_SFo2i%copy_only) then
               call seq_map_set_type(mapper_SFo2i, mboxid, 1) ! type is cells
            endif

         endif
#endif
       endif

       if (glc_c2_ice) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rg2i'
          end if
          call seq_map_init_rcfile(mapper_Rg2i, glc(1), ice(1), &
               'seq_maps.rc','glc2ice_rmapname:','glc2ice_rmaptype:',samegrid_ig, &
               'mapper_Rg2i initialization', esmf_map_flag)
       endif

       if (glcshelf_c2_ice) then
         if (iamroot_CPLID) then
               write(logunit,*) ' '
               write(logunit,F00) 'Initializing mapper_Sg2i'
         end if
         call seq_map_init_rcfile(mapper_Sg2i, glc(1), ice(1), &
              'seq_maps.rc','glc2ice_smapname:','glc2ice_smaptype:',samegrid_ig, &
              'mapper_Sg2i initialization', esmf_map_flag)

         if (iamroot_CPLID) then
               write(logunit,*) ' '
               write(logunit,F00) 'Initializing mapper_Fg2i'
         end if
         call seq_map_init_rcfile(mapper_Fg2i, glc(1), ice(1), &
              'seq_maps.rc','glc2ice_fmapname:','glc2ice_fmaptype:',samegrid_ig, &
              'mapper_Fg2i initialization', esmf_map_flag)
       endif

       if (rof_c2_ice) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rr2i'
          end if
          call seq_map_init_rcfile(mapper_Rr2i, rof(1), ice(1), &
               'seq_maps.rc','rof2ice_rmapname:','rof2ice_rmaptype:',samegrid_ro, &
               'mapper_Rr2i initialization', esmf_map_flag)
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_ice_init

  !================================================================================================

  subroutine prep_ice_mrg(infodata, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Prepare run phase, including running the merge
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eoi, eai, egi, eii, eri
    real(r8)                 :: flux_epbalfact ! adjusted precip factor
    type(mct_avect), pointer :: x2i_ix
    character(*), parameter  :: subname = '(prep_ice_mrg)'
    !---------------------------------------------------------------

    call seq_infodata_GetData(infodata, &
         flux_epbalfact=flux_epbalfact)

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do eii = 1,num_inst_ice
       ! Use fortran mod to address ensembles in merge
       eai = mod((eii-1),num_inst_atm) + 1
       eoi = mod((eii-1),num_inst_ocn) + 1
       eri = mod((eii-1),num_inst_rof) + 1
       egi = mod((eii-1),num_inst_glc) + 1

       ! Apply correction to precipitation of requested driver namelist
       x2i_ix   => component_get_x2c_cx(ice(eii))  ! This is actually modifying x2i_ix
       call prep_ice_merge(flux_epbalfact, a2x_ix(eai), o2x_ix(eoi), r2x_ix(eri), g2x_ix(egi), &
            x2i_ix)
    enddo
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_ice_mrg

  !================================================================================================

  subroutine prep_ice_merge(flux_epbalfact, a2x_i, o2x_i, r2x_i, g2x_i, x2i_i )

    !-----------------------------------------------------------------------
    !
    ! Arguments
    real(r8)        , intent(inout) :: flux_epbalfact
    type(mct_aVect) , intent(in)    :: a2x_i
    type(mct_aVect) , intent(in)    :: o2x_i
    type(mct_aVect) , intent(in)    :: r2x_i
    type(mct_aVect) , intent(in)    :: g2x_i
    type(mct_aVect) , intent(inout) :: x2i_i
    !
    ! Local variables
    integer       :: i,i1,o1,lsize
    integer       :: niflds
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_g2x_Figg_rofi
    integer, save :: index_r2x_Firr_rofi
    integer, save :: index_x2i_Faxa_rain
    integer, save :: index_x2i_Faxa_snow
    integer, save :: index_x2i_Fixx_rofi
    !wiso fields:
    integer, save :: index_a2x_Faxa_rainc_16O
    integer, save :: index_a2x_Faxa_rainl_16O
    integer, save :: index_a2x_Faxa_snowc_16O
    integer, save :: index_a2x_Faxa_snowl_16O
    integer, save :: index_x2i_Faxa_rain_16O
    integer, save :: index_x2i_Faxa_snow_16O
    integer, save :: index_a2x_Faxa_rainc_18O
    integer, save :: index_a2x_Faxa_rainl_18O
    integer, save :: index_a2x_Faxa_snowc_18O
    integer, save :: index_a2x_Faxa_snowl_18O
    integer, save :: index_x2i_Faxa_rain_18O
    integer, save :: index_x2i_Faxa_snow_18O
    integer, save :: index_a2x_Faxa_rainc_HDO
    integer, save :: index_a2x_Faxa_rainl_HDO
    integer, save :: index_a2x_Faxa_snowc_HDO
    integer, save :: index_a2x_Faxa_snowl_HDO
    integer, save :: index_x2i_Faxa_rain_HDO
    integer, save :: index_x2i_Faxa_snow_HDO
    logical, save :: first_time = .true.
    logical       :: iamroot
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: o2x_sharedindices
    type(mct_aVect_sharedindices),save :: a2x_sharedindices
    type(mct_aVect_sharedindices),save :: g2x_sharedindices
    character(*), parameter   :: subname = '(prep_ice_merge) '
    !-----------------------------------------------------------------------

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2i_i)

    if (first_time) then
       niflds = mct_aVect_nRattr(x2i_i)

       allocate(mrgstr(niflds))
       index_a2x_Faxa_snowc = mct_aVect_indexRA(a2x_i,'Faxa_snowc')
       index_a2x_Faxa_snowl = mct_aVect_indexRA(a2x_i,'Faxa_snowl')
       index_a2x_Faxa_rainc = mct_aVect_indexRA(a2x_i,'Faxa_rainc')
       index_a2x_Faxa_rainl = mct_aVect_indexRA(a2x_i,'Faxa_rainl')
       index_g2x_Figg_rofi  = mct_aVect_indexRA(g2x_i,'Figg_rofi')
       index_r2x_Firr_rofi  = mct_aVect_indexRA(r2x_i,'Firr_rofi')
       index_x2i_Faxa_rain  = mct_aVect_indexRA(x2i_i,'Faxa_rain' )
       index_x2i_Faxa_snow  = mct_aVect_indexRA(x2i_i,'Faxa_snow' )
       index_x2i_Fixx_rofi  = mct_aVect_indexRA(x2i_i,'Fixx_rofi')

       ! Water isotope fields
       index_a2x_Faxa_snowc_16O = mct_aVect_indexRA(a2x_i,'Faxa_snowc_16O', perrWith='quiet')
       index_a2x_Faxa_snowl_16O = mct_aVect_indexRA(a2x_i,'Faxa_snowl_16O', perrWith='quiet')
       index_a2x_Faxa_rainc_16O = mct_aVect_indexRA(a2x_i,'Faxa_rainc_16O', perrWith='quiet')
       index_a2x_Faxa_rainl_16O = mct_aVect_indexRA(a2x_i,'Faxa_rainl_16O', perrWith='quiet')
       index_x2i_Faxa_rain_16O  = mct_aVect_indexRA(x2i_i,'Faxa_rain_16O',  perrWith='quiet' )
       index_x2i_Faxa_snow_16O  = mct_aVect_indexRA(x2i_i,'Faxa_snow_16O',  perrWith='quiet' )

       index_a2x_Faxa_snowc_18O = mct_aVect_indexRA(a2x_i,'Faxa_snowc_18O', perrWith='quiet')
       index_a2x_Faxa_snowl_18O = mct_aVect_indexRA(a2x_i,'Faxa_snowl_18O', perrWith='quiet')
       index_a2x_Faxa_rainc_18O = mct_aVect_indexRA(a2x_i,'Faxa_rainc_18O', perrWith='quiet')
       index_a2x_Faxa_rainl_18O = mct_aVect_indexRA(a2x_i,'Faxa_rainl_18O', perrWith='quiet')
       index_x2i_Faxa_rain_18O  = mct_aVect_indexRA(x2i_i,'Faxa_rain_18O',  perrWith='quiet' )
       index_x2i_Faxa_snow_18O  = mct_aVect_indexRA(x2i_i,'Faxa_snow_18O',  perrWith='quiet' )

       index_a2x_Faxa_snowc_HDO = mct_aVect_indexRA(a2x_i,'Faxa_snowc_HDO', perrWith='quiet')
       index_a2x_Faxa_snowl_HDO = mct_aVect_indexRA(a2x_i,'Faxa_snowl_HDO', perrWith='quiet')
       index_a2x_Faxa_rainc_HDO = mct_aVect_indexRA(a2x_i,'Faxa_rainc_HDO', perrWith='quiet')
       index_a2x_Faxa_rainl_HDO = mct_aVect_indexRA(a2x_i,'Faxa_rainl_HDO', perrWith='quiet')
       index_x2i_Faxa_rain_HDO  = mct_aVect_indexRA(x2i_i,'Faxa_rain_HDO',  perrWith='quiet' )
       index_x2i_Faxa_snow_HDO  = mct_aVect_indexRA(x2i_i,'Faxa_snow_HDO',  perrWith='quiet' )

       do i = 1,niflds
          field = mct_aVect_getRList2c(i, x2i_i)
          mrgstr(i) = subname//'x2i%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(o2x_i, x2i_i, o2x_SharedIndices)
       call mct_aVect_setSharedIndices(a2x_i, x2i_i, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(g2x_i, x2i_i, g2x_SharedIndices)

       !--- document copy operations ---
       do i=1,o2x_SharedIndices%shared_real%num_indices
          i1=o2x_SharedIndices%shared_real%aVindices1(i)
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, o2x_i)
          mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field)
       enddo
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, a2x_i)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
       enddo
       do i=1,g2x_SharedIndices%shared_real%num_indices
          i1=g2x_SharedIndices%shared_real%aVindices1(i)
          o1=g2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, g2x_i)
          mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field)
       enddo

       !--- document manual merges ---
       mrgstr(index_x2i_Faxa_rain) = trim(mrgstr(index_x2i_Faxa_rain))//' = '// &
            '(a2x%Faxa_rainc + a2x%Faxa_rainl)*flux_epbalfact'
       mrgstr(index_x2i_Faxa_snow) = trim(mrgstr(index_x2i_Faxa_snow))//' = '// &
            '(a2x%Faxa_snowc + a2x%Faxa_snowl)*flux_epbalfact'
       mrgstr(index_x2i_Fixx_rofi) = trim(mrgstr(index_x2i_Fixx_rofi))//' = '// &
            '(g2x%Figg_rofi + r2x%Firr_rofi)*flux_epbalfact'

       !--- water isotope document manual merges ---
       if ( index_x2i_Faxa_rain_16O /= 0 ) then
          mrgstr(index_x2i_Faxa_rain_16O) = trim(mrgstr(index_x2i_Faxa_rain_16O))//' = '// &
               '(a2x%Faxa_rainc_16O + a2x%Faxa_rainl_16O)*flux_epbalfact'
          mrgstr(index_x2i_Faxa_snow_16O) = trim(mrgstr(index_x2i_Faxa_snow_16O))//' = '// &
               '(a2x%Faxa_snowc_16O + a2x%Faxa_snowl_16O)*flux_epbalfact'
       end if
       if ( index_x2i_Faxa_rain_18O /= 0 ) then
          mrgstr(index_x2i_Faxa_rain_18O) = trim(mrgstr(index_x2i_Faxa_rain_18O))//' = '// &
               '(a2x%Faxa_rainc_18O + a2x%Faxa_rainl_18O)*flux_epbalfact'
          mrgstr(index_x2i_Faxa_snow_18O) = trim(mrgstr(index_x2i_Faxa_snow_18O))//' = '// &
               '(a2x%Faxa_snowc_18O + a2x%Faxa_snowl_18O)*flux_epbalfact'
       end if
       if ( index_x2i_Faxa_rain_HDO /= 0 ) then
          mrgstr(index_x2i_Faxa_rain_HDO) = trim(mrgstr(index_x2i_Faxa_rain_HDO))//' = '// &
               '(a2x%Faxa_rainc_HDO + a2x%Faxa_rainl_HDO)*flux_epbalfact'
          mrgstr(index_x2i_Faxa_snow_HDO) = trim(mrgstr(index_x2i_Faxa_snow_HDO))//' = '// &
               '(a2x%Faxa_snowc_HDO + a2x%Faxa_snowl_HDO)*flux_epbalfact'
       end if

    endif

    !    call mct_aVect_copy(aVin=o2x_i, aVout=x2i_i, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=a2x_i, aVout=x2i_i, vector=mct_usevector)
    !    call mct_aVect_copy(aVin=g2x_i, aVout=x2i_i, vector=mct_usevector)
    call mct_aVect_copy(aVin=o2x_i, aVout=x2i_i, vector=mct_usevector, sharedIndices=o2x_SharedIndices)
    call mct_aVect_copy(aVin=a2x_i, aVout=x2i_i, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    call mct_aVect_copy(aVin=g2x_i, aVout=x2i_i, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

    ! Merge total snow and precip for ice input
    ! Scale total precip and runoff by flux_epbalfact

    do i = 1,lsize
       x2i_i%rAttr(index_x2i_Faxa_rain,i) = a2x_i%rAttr(index_a2x_Faxa_rainc,i) + &
            a2x_i%rAttr(index_a2x_Faxa_rainl,i)
       x2i_i%rAttr(index_x2i_Faxa_snow,i) = a2x_i%rAttr(index_a2x_Faxa_snowc,i) + &
            a2x_i%rAttr(index_a2x_Faxa_snowl,i)
       x2i_i%rAttr(index_x2i_Fixx_rofi,i) = g2x_i%rAttr(index_g2x_Figg_rofi,i) + &
            r2x_i%rAttr(index_r2x_Firr_rofi,i)

       x2i_i%rAttr(index_x2i_Faxa_rain,i) = x2i_i%rAttr(index_x2i_Faxa_rain,i) * flux_epbalfact
       x2i_i%rAttr(index_x2i_Faxa_snow,i) = x2i_i%rAttr(index_x2i_Faxa_snow,i) * flux_epbalfact
       x2i_i%rAttr(index_x2i_Fixx_rofi,i) = x2i_i%rAttr(index_x2i_Fixx_rofi,i) * flux_epbalfact

       ! For water isotopes
       if ( index_x2i_Faxa_rain_16O /= 0 ) then
          x2i_i%rAttr(index_x2i_Faxa_rain_16O,i) = a2x_i%rAttr(index_a2x_Faxa_rainc_16O,i) + &
               a2x_i%rAttr(index_a2x_Faxa_rainl_16O,i)
          x2i_i%rAttr(index_x2i_Faxa_snow_16O,i) = a2x_i%rAttr(index_a2x_Faxa_snowc_16O,i) + &
               a2x_i%rAttr(index_a2x_Faxa_snowl_16O,i)

          x2i_i%rAttr(index_x2i_Faxa_rain_16O,i) = x2i_i%rAttr(index_x2i_Faxa_rain_16O,i) * flux_epbalfact
          x2i_i%rAttr(index_x2i_Faxa_snow_16O,i) = x2i_i%rAttr(index_x2i_Faxa_snow_16O,i) * flux_epbalfact
       end if
       if ( index_x2i_Faxa_rain_18O /= 0 ) then
          x2i_i%rAttr(index_x2i_Faxa_rain_18O,i) = a2x_i%rAttr(index_a2x_Faxa_rainc_18O,i) + &
               a2x_i%rAttr(index_a2x_Faxa_rainl_18O,i)
          x2i_i%rAttr(index_x2i_Faxa_snow_18O,i) = a2x_i%rAttr(index_a2x_Faxa_snowc_18O,i) + &
               a2x_i%rAttr(index_a2x_Faxa_snowl_18O,i)

          x2i_i%rAttr(index_x2i_Faxa_rain_18O,i) = x2i_i%rAttr(index_x2i_Faxa_rain_18O,i) * flux_epbalfact
          x2i_i%rAttr(index_x2i_Faxa_snow_18O,i) = x2i_i%rAttr(index_x2i_Faxa_snow_18O,i) * flux_epbalfact
       end if
       if ( index_x2i_Faxa_rain_HDO /= 0 ) then
          x2i_i%rAttr(index_x2i_Faxa_rain_HDO,i) = a2x_i%rAttr(index_a2x_Faxa_rainc_HDO,i) + &
               a2x_i%rAttr(index_a2x_Faxa_rainl_HDO,i)
          x2i_i%rAttr(index_x2i_Faxa_snow_HDO,i) = a2x_i%rAttr(index_a2x_Faxa_snowc_HDO,i) + &
               a2x_i%rAttr(index_a2x_Faxa_snowl_HDO,i)

          x2i_i%rAttr(index_x2i_Faxa_rain_HDO,i) = x2i_i%rAttr(index_x2i_Faxa_rain_HDO,i) * flux_epbalfact
          x2i_i%rAttr(index_x2i_Faxa_snow_HDO,i) = x2i_i%rAttr(index_x2i_Faxa_snow_HDO,i) * flux_epbalfact
       end if

    end do

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,niflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif

    first_time = .false.

  end subroutine prep_ice_merge

  subroutine prep_ice_mrg_moab(infodata, rof_c2_ice)
    use iMOAB , only : iMOAB_GetDoubleTagStorage, &
    iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh, iMOAB_GetMeshInfo

    !-----------------------------------------------------------------------
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical,                   intent(in)    :: rof_c2_ice ! .true.  => rof to ice coupling on
    !
    ! Local variables
    real(r8)        :: flux_epbalfact
    integer       :: i,i1,o1, arrsize
    integer,save  :: niflds,naflds,nrflds
    integer       :: lsize
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_g2x_Figg_rofi
    integer, save :: index_r2x_Firr_rofi
    integer, save :: index_x2i_Faxa_rain
    integer, save :: index_x2i_Faxa_snow
    integer, save :: index_x2i_Fixx_rofi
    !wiso fields:
    integer, save :: index_a2x_Faxa_rainc_16O
    integer, save :: index_a2x_Faxa_rainl_16O
    integer, save :: index_a2x_Faxa_snowc_16O
    integer, save :: index_a2x_Faxa_snowl_16O
    integer, save :: index_x2i_Faxa_rain_16O
    integer, save :: index_x2i_Faxa_snow_16O
    integer, save :: index_a2x_Faxa_rainc_18O
    integer, save :: index_a2x_Faxa_rainl_18O
    integer, save :: index_a2x_Faxa_snowc_18O
    integer, save :: index_a2x_Faxa_snowl_18O
    integer, save :: index_x2i_Faxa_rain_18O
    integer, save :: index_x2i_Faxa_snow_18O
    integer, save :: index_a2x_Faxa_rainc_HDO
    integer, save :: index_a2x_Faxa_rainl_HDO
    integer, save :: index_a2x_Faxa_snowc_HDO
    integer, save :: index_a2x_Faxa_snowl_HDO
    integer, save :: index_x2i_Faxa_rain_HDO
    integer, save :: index_x2i_Faxa_snow_HDO
    logical, save :: first_time = .true.
    logical       :: iamroot

    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info
    character(CXX) ::tagname, mct_field
    character(CL),allocatable :: mrgstr(:)   ! temporary string

    character(CL) :: field   ! string converted to char

    type(mct_avect) , pointer, save   :: a2x_i
    type(mct_avect) , pointer, save   :: x2i_i
    type(mct_avect) , pointer, save   :: r2x_i
    type(mct_avect) , pointer, save   :: o2x_i
    type(mct_avect) , pointer, save   :: g2x_i
    type(mct_aVect_sharedindices),save :: o2x_SharedIndices
    type(mct_aVect_sharedindices),save :: a2x_SharedIndices
    type(mct_aVect_sharedindices),save :: g2x_SharedIndices

    integer ent_type, ierr,n
#ifdef MOABDEBUG
    character*32             :: outfile, wopts, lnum
#endif
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
#endif


    character(*), parameter   :: subname = '(prep_ice_mrg_moab) '
    !-----------------------------------------------------------------------
    call seq_infodata_GetData(infodata, &
         flux_epbalfact=flux_epbalfact)

    call seq_comm_getdata(CPLID, iamroot=iamroot)

 ! find out the number of local elements in moab mesh seaice instance on coupler
    ierr  = iMOAB_GetMeshInfo ( mbixid, nvert, nvise, nbl, nsurf, nvisBC );
    if (ierr .ne. 0) then
         write(logunit,*) subname,' error in getting info '
         call shr_sys_abort(subname//' error in getting info ')
    endif
    lsize = nvise(1) ! number of active cells

    if (first_time) then
      a2x_i => a2x_ix(1)
      r2x_i => r2x_ix(1)
      o2x_i => o2x_ix(1)
      g2x_i => g2x_ix(1)
      x2i_i => component_get_x2c_cx(ice(1))

      niflds = mct_aVect_nRattr(x2i_i)
      naflds = mct_aVect_nRattr(a2x_i)
      nrflds = mct_aVect_nRattr(r2x_i)


       allocate(x2i_im(lsize,niflds))
       allocate(a2x_im(lsize,naflds))
       allocate(r2x_im(lsize,nrflds))
       r2x_im = 0._r8 ! should we zero out all of them ?

       allocate(mrgstr(niflds))
       index_a2x_Faxa_snowc = mct_aVect_indexRA(a2x_i,'Faxa_snowc')
       index_a2x_Faxa_snowl = mct_aVect_indexRA(a2x_i,'Faxa_snowl')
       index_a2x_Faxa_rainc = mct_aVect_indexRA(a2x_i,'Faxa_rainc')
       index_a2x_Faxa_rainl = mct_aVect_indexRA(a2x_i,'Faxa_rainl')
       index_g2x_Figg_rofi  = mct_aVect_indexRA(g2x_i,'Figg_rofi')
       index_r2x_Firr_rofi  = mct_aVect_indexRA(r2x_i,'Firr_rofi')
       index_x2i_Faxa_rain  = mct_aVect_indexRA(x2i_i,'Faxa_rain' )
       index_x2i_Faxa_snow  = mct_aVect_indexRA(x2i_i,'Faxa_snow' )
       index_x2i_Fixx_rofi  = mct_aVect_indexRA(x2i_i,'Fixx_rofi')

       ! Water isotope fields
       index_a2x_Faxa_snowc_16O = mct_aVect_indexRA(a2x_i,'Faxa_snowc_16O', perrWith='quiet')
       index_a2x_Faxa_snowl_16O = mct_aVect_indexRA(a2x_i,'Faxa_snowl_16O', perrWith='quiet')
       index_a2x_Faxa_rainc_16O = mct_aVect_indexRA(a2x_i,'Faxa_rainc_16O', perrWith='quiet')
       index_a2x_Faxa_rainl_16O = mct_aVect_indexRA(a2x_i,'Faxa_rainl_16O', perrWith='quiet')
       index_x2i_Faxa_rain_16O  = mct_aVect_indexRA(x2i_i,'Faxa_rain_16O',  perrWith='quiet' )
       index_x2i_Faxa_snow_16O  = mct_aVect_indexRA(x2i_i,'Faxa_snow_16O',  perrWith='quiet' )

       index_a2x_Faxa_snowc_18O = mct_aVect_indexRA(a2x_i,'Faxa_snowc_18O', perrWith='quiet')
       index_a2x_Faxa_snowl_18O = mct_aVect_indexRA(a2x_i,'Faxa_snowl_18O', perrWith='quiet')
       index_a2x_Faxa_rainc_18O = mct_aVect_indexRA(a2x_i,'Faxa_rainc_18O', perrWith='quiet')
       index_a2x_Faxa_rainl_18O = mct_aVect_indexRA(a2x_i,'Faxa_rainl_18O', perrWith='quiet')
       index_x2i_Faxa_rain_18O  = mct_aVect_indexRA(x2i_i,'Faxa_rain_18O',  perrWith='quiet' )
       index_x2i_Faxa_snow_18O  = mct_aVect_indexRA(x2i_i,'Faxa_snow_18O',  perrWith='quiet' )

       index_a2x_Faxa_snowc_HDO = mct_aVect_indexRA(a2x_i,'Faxa_snowc_HDO', perrWith='quiet')
       index_a2x_Faxa_snowl_HDO = mct_aVect_indexRA(a2x_i,'Faxa_snowl_HDO', perrWith='quiet')
       index_a2x_Faxa_rainc_HDO = mct_aVect_indexRA(a2x_i,'Faxa_rainc_HDO', perrWith='quiet')
       index_a2x_Faxa_rainl_HDO = mct_aVect_indexRA(a2x_i,'Faxa_rainl_HDO', perrWith='quiet')
       index_x2i_Faxa_rain_HDO  = mct_aVect_indexRA(x2i_i,'Faxa_rain_HDO',  perrWith='quiet' )
       index_x2i_Faxa_snow_HDO  = mct_aVect_indexRA(x2i_i,'Faxa_snow_HDO',  perrWith='quiet' )

       do i = 1,niflds
          field = mct_aVect_getRList2c(i, x2i_i)
          mrgstr(i) = subname//'x2i%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(o2x_i, x2i_i, o2x_SharedIndices)
       call mct_aVect_setSharedIndices(a2x_i, x2i_i, a2x_SharedIndices)
       call mct_aVect_setSharedIndices(g2x_i, x2i_i, g2x_SharedIndices)

       !--- document copy operations ---
       do i=1,o2x_SharedIndices%shared_real%num_indices
          i1=o2x_SharedIndices%shared_real%aVindices1(i)
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, o2x_i)
          mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field)
#ifdef MOABDEBUG
          write(lnum, "(I3, A18, I3)" )i1, ' in o2x_ix, x2i_ix ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,a2x_SharedIndices%shared_real%num_indices
          i1=a2x_SharedIndices%shared_real%aVindices1(i)
          o1=a2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, a2x_i)
          mrgstr(o1) = trim(mrgstr(o1))//' = a2x%'//trim(field)
#ifdef MOABDEBUG
          write(lnum, "(I3, A18, I3)" )i1, ' in a2x_ix, x2i_ix', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo
       do i=1,g2x_SharedIndices%shared_real%num_indices
          i1=g2x_SharedIndices%shared_real%aVindices1(i)
          o1=g2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, g2x_i)
          mrgstr(o1) = trim(mrgstr(o1))//' = g2x%'//trim(field)
#ifdef MOABDEBUG
          write(lnum, "(I3, A18, I3)" )i1, ' in g2x_ix, x2i_ix ', o1
          mrgstr(o1) = trim(mrgstr(o1))//trim(lnum)
#endif
       enddo

       !--- document manual merges ---
       mrgstr(index_x2i_Faxa_rain) = trim(mrgstr(index_x2i_Faxa_rain))//' = '// &
            '(a2x%Faxa_rainc + a2x%Faxa_rainl)*flux_epbalfact'
       mrgstr(index_x2i_Faxa_snow) = trim(mrgstr(index_x2i_Faxa_snow))//' = '// &
            '(a2x%Faxa_snowc + a2x%Faxa_snowl)*flux_epbalfact'
       mrgstr(index_x2i_Fixx_rofi) = trim(mrgstr(index_x2i_Fixx_rofi))//' = '// &
            '(g2x%Figg_rofi + r2x%Firr_rofi)*flux_epbalfact'

       !--- water isotope document manual merges ---
       if ( index_x2i_Faxa_rain_16O /= 0 ) then
          mrgstr(index_x2i_Faxa_rain_16O) = trim(mrgstr(index_x2i_Faxa_rain_16O))//' = '// &
               '(a2x%Faxa_rainc_16O + a2x%Faxa_rainl_16O)*flux_epbalfact'
          mrgstr(index_x2i_Faxa_snow_16O) = trim(mrgstr(index_x2i_Faxa_snow_16O))//' = '// &
               '(a2x%Faxa_snowc_16O + a2x%Faxa_snowl_16O)*flux_epbalfact'
       end if
       if ( index_x2i_Faxa_rain_18O /= 0 ) then
          mrgstr(index_x2i_Faxa_rain_18O) = trim(mrgstr(index_x2i_Faxa_rain_18O))//' = '// &
               '(a2x%Faxa_rainc_18O + a2x%Faxa_rainl_18O)*flux_epbalfact'
          mrgstr(index_x2i_Faxa_snow_18O) = trim(mrgstr(index_x2i_Faxa_snow_18O))//' = '// &
               '(a2x%Faxa_snowc_18O + a2x%Faxa_snowl_18O)*flux_epbalfact'
       end if
       if ( index_x2i_Faxa_rain_HDO /= 0 ) then
          mrgstr(index_x2i_Faxa_rain_HDO) = trim(mrgstr(index_x2i_Faxa_rain_HDO))//' = '// &
               '(a2x%Faxa_rainc_HDO + a2x%Faxa_rainl_HDO)*flux_epbalfact'
          mrgstr(index_x2i_Faxa_snow_HDO) = trim(mrgstr(index_x2i_Faxa_snow_HDO))//' = '// &
               '(a2x%Faxa_snowc_HDO + a2x%Faxa_snowl_HDO)*flux_epbalfact'
       end if

    endif

    !call mct_aVect_copy(aVin=o2x_i, aVout=x2i_i, vector=mct_usevector, sharedIndices=o2x_SharedIndices)
    !call mct_aVect_copy(aVin=a2x_i, aVout=x2i_i, vector=mct_usevector, sharedIndices=a2x_SharedIndices)
    !call mct_aVect_copy(aVin=g2x_i, aVout=x2i_i, vector=mct_usevector, sharedIndices=g2x_SharedIndices)

    ! get data that needs custom merges
    ent_type = 1

   ! git the x2i_im field which has been mostly filled out by mapping calls.
    tagname = trim(seq_flds_x2i_fields)//C_NULL_CHAR
    arrsize = niflds * lsize
    ierr = iMOAB_GetDoubleTagStorage ( mbixid, tagname, arrsize , ent_type, x2i_im)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting x2i_im array ')
    endif

! get the a2x data that was mapped to i
    tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
    arrsize = naflds * lsize !        allocate (a2x_om (lsize, naflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbixid, tagname, arrsize , ent_type, a2x_im)
    if (ierr .ne. 0) then
      write(logunit, *) 'MOAB error ', ierr
      call shr_sys_abort(subname//' error in getting a2x_im array ')
    endif

! get the r2x data that was mapped to i
    tagname = trim(seq_flds_r2x_fields)//C_NULL_CHAR
    arrsize = nrflds * lsize !        allocate (a2x_om (lsize, naflds))
    ierr = iMOAB_GetDoubleTagStorage ( mbixid, tagname, arrsize , ent_type, r2x_im)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in getting r2x_im array ')
    endif

    ! Merge total snow and precip for ice input
    ! Scale total precip and runoff by flux_epbalfact


    do n = 1,lsize
       x2i_im(n,index_x2i_Faxa_rain) = a2x_im(n,index_a2x_Faxa_rainc) + a2x_im(n,index_a2x_Faxa_rainl)
       x2i_im(n,index_x2i_Faxa_snow) = a2x_im(n,index_a2x_Faxa_snowc) + a2x_im(n,index_a2x_Faxa_snowl)

! no glacier yet
!       x2i_im(n,index_x2i_Fixx_rofi) = g2x_im(n,index_g2x_Figg_rofi) + &
!            r2x_im(n,index_r2x_Firr_rofi)
        if(rof_c2_ice) then
            x2i_im(n,index_x2i_Fixx_rofi) = r2x_im(n,index_r2x_Firr_rofi)
        else
            x2i_im(n,index_x2i_Fixx_rofi) = 0._r8
        endif

       x2i_im(n,index_x2i_Faxa_rain) = x2i_im(n,index_x2i_Faxa_rain) * flux_epbalfact
       x2i_im(n,index_x2i_Faxa_snow) = x2i_im(n,index_x2i_Faxa_snow) * flux_epbalfact
       x2i_im(n,index_x2i_Fixx_rofi) = x2i_im(n,index_x2i_Fixx_rofi) * flux_epbalfact

       ! For water isotopes
       if ( index_x2i_Faxa_rain_16O /= 0 ) then
          x2i_im(n,index_x2i_Faxa_rain_16O) = a2x_im(n,index_a2x_Faxa_rainc_16O) + &
               a2x_im(n,index_a2x_Faxa_rainl_16O)
          x2i_im(n,index_x2i_Faxa_snow_16O) = a2x_im(n,index_a2x_Faxa_snowc_16O) + &
               a2x_im(n,index_a2x_Faxa_snowl_16O)

          x2i_im(n,index_x2i_Faxa_rain_16O) = x2i_im(n,index_x2i_Faxa_rain_16O) * flux_epbalfact
          x2i_im(n,index_x2i_Faxa_snow_16O) = x2i_im(n,index_x2i_Faxa_snow_16O) * flux_epbalfact
       end if
       if ( index_x2i_Faxa_rain_18O /= 0 ) then
          x2i_im(n,index_x2i_Faxa_rain_18O) = a2x_im(n,index_a2x_Faxa_rainc_18O) + &
               a2x_im(n,index_a2x_Faxa_rainl_18O)
          x2i_im(n,index_x2i_Faxa_snow_18O) = a2x_im(n,index_a2x_Faxa_snowc_18O) + &
               a2x_im(n,index_a2x_Faxa_snowl_18O)

          x2i_im(n,index_x2i_Faxa_rain_18O) = x2i_im(n,index_x2i_Faxa_rain_18O) * flux_epbalfact
          x2i_im(n,index_x2i_Faxa_snow_18O) = x2i_im(n,index_x2i_Faxa_snow_18O) * flux_epbalfact
       end if
       if ( index_x2i_Faxa_rain_HDO /= 0 ) then
          x2i_im(n,index_x2i_Faxa_rain_HDO) = a2x_im(n,index_a2x_Faxa_rainc_HDO) + &
               a2x_im(n,index_a2x_Faxa_rainl_HDO)
          x2i_im(n,index_x2i_Faxa_snow_HDO) = a2x_im(n,index_a2x_Faxa_snowc_HDO) + &
               a2x_im(n,index_a2x_Faxa_snowl_HDO)

          x2i_im(n,index_x2i_Faxa_rain_HDO) = x2i_im(n,index_x2i_Faxa_rain_HDO) * flux_epbalfact
          x2i_im(n,index_x2i_Faxa_snow_HDO) = x2i_im(n,index_x2i_Faxa_snow_HDO) * flux_epbalfact
       end if

    end do
    tagname = trim(seq_flds_x2i_fields)//C_NULL_CHAR
    arrsize = niflds * lsize
    ierr = iMOAB_SetDoubleTagStorage ( mbixid, tagname, arrsize , ent_type, x2i_im)
    if (ierr .ne. 0) then
      call shr_sys_abort(subname//' error in setting x2i_im array ')
    endif

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary MOAB:'
          do i = 1,niflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif



#ifdef MOABCOMP
 !compare_mct_av_moab_tag(comp, attrVect, field, imoabApp, tag_name, ent_type, difference)
    x2i_i => component_get_x2c_cx(ice(1))
    ! loop over all fields in seq_flds_x2i_fields
    call mct_list_init(temp_list ,seq_flds_x2i_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! cell for ice/ocean
    if (iamroot) print *, num_moab_exports, trim(seq_flds_x2i_fields)
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call compare_mct_av_moab_tag(ice(1), x2i_i, mct_field,  mbixid, tagname, ent_type, difference, first_time)
    enddo
    call mct_list_clean(temp_list)
#endif
    first_time = .false.
#ifdef MOABDEBUG
    if (mbixid .ge. 0 ) then !  we are on coupler pes, for sure
       write(lnum,"(I0.2)")num_moab_exports
       outfile = 'IceCplAftMm'//trim(lnum)//'.h5m'//C_NULL_CHAR
       wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
       ierr = iMOAB_WriteMesh(mbixid, trim(outfile), trim(wopts))
    endif
#endif


  end subroutine prep_ice_mrg_moab

  subroutine prep_ice_calc_a2x_ix(a2x_ox, timer)
    use iMOAB , only : iMOAB_WriteMesh
    !---------------------------------------------------------------
    ! Description
    ! Create a2x_ix (note that a2x_ix is a local module variable)
    !
    ! Arguments
    type(mct_aVect) , intent(in) :: a2x_ox(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    character(*), parameter :: subname = '(prep_ice_calc_a2x_ix)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       call seq_map_map(mapper_SFo2i, a2x_ox(eai), a2x_ix(eai), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ice_calc_a2x_ix

  !================================================================================================

  subroutine prep_ice_calc_o2x_ix(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_ix (note that o2x_ix is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi
    type(mct_aVect) , pointer :: o2x_ox
    character(*), parameter :: subname = '(prep_ice_calc_o2x_ix)'

    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       o2x_ox => component_get_c2x_cx(ocn(eoi))
       call seq_map_map(mapper_SFo2i, o2x_ox, o2x_ix(eoi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ice_calc_o2x_ix

  !================================================================================================

  subroutine prep_ice_calc_r2x_ix(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create r2x_ix (note that r2x_ix is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri
    type(mct_aVect), pointer :: r2x_rx
    character(*), parameter :: subname = '(prep_ice_calc_r2x_ix)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       r2x_rx => component_get_c2x_cx(rof(eri))

       call seq_map_map(mapper_Rr2i, r2x_rx, r2x_ix(eri), &
            fldlist='Firr_rofi', norm=.false.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ice_calc_r2x_ix

  !================================================================================================

  subroutine prep_ice_calc_g2x_ix(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_ix (note that g2x_ix is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_aVect), pointer :: g2x_gx
    character(*), parameter :: subname = '(prep_ice_calc_g2x_ix)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       g2x_gx => component_get_c2x_cx(glc(egi))
       call seq_map_map(mapper_Rg2i, g2x_gx, g2x_ix(egi), &
                        fldlist='Fixx_rofi', norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ice_calc_g2x_ix

  !================================================================================================

  subroutine prep_ice_shelf_calc_g2x_ix(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_ix (note that g2x_ix is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_aVect), pointer :: g2x_gx
    character(*), parameter :: subname = '(prep_ice_calc_g2x_ix)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_rof
       g2x_gx => component_get_c2x_cx(glc(egi))
       call seq_map_map(mapper_Sg2i, g2x_gx, g2x_ix(egi), &
                        fldlist='Sg_icemask_coupled_fluxes', norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ice_shelf_calc_g2x_ix

  !================================================================================================

  function prep_ice_get_a2x_ix()
    type(mct_aVect), pointer :: prep_ice_get_a2x_ix(:)
    prep_ice_get_a2x_ix => a2x_ix(:)
  end function prep_ice_get_a2x_ix

  function prep_ice_get_o2x_ix()
    type(mct_aVect), pointer :: prep_ice_get_o2x_ix(:)
    prep_ice_get_o2x_ix => o2x_ix(:)
  end function prep_ice_get_o2x_ix

  function prep_ice_get_g2x_ix()
    type(mct_aVect), pointer :: prep_ice_get_g2x_ix(:)
    prep_ice_get_g2x_ix => g2x_ix(:)
  end function prep_ice_get_g2x_ix

  function prep_ice_get_r2x_ix()
    type(mct_aVect), pointer :: prep_ice_get_r2x_ix(:)
    prep_ice_get_r2x_ix => r2x_ix(:)
  end function prep_ice_get_r2x_ix

  function prep_ice_get_mapper_SFo2i()
    type(seq_map), pointer :: prep_ice_get_mapper_SFo2i
    prep_ice_get_mapper_SFo2i => mapper_SFo2i
  end function prep_ice_get_mapper_SFo2i

  function prep_ice_get_mapper_Rg2i()
    type(seq_map), pointer :: prep_ice_get_mapper_Rg2i
    prep_ice_get_mapper_Rg2i => mapper_Rg2i
  end function prep_ice_get_mapper_Rg2i

  function prep_ice_get_mapper_Sg2i()
    type(seq_map), pointer :: prep_ice_get_mapper_Sg2i
    prep_ice_get_mapper_Sg2i => mapper_Sg2i
  end function prep_ice_get_mapper_Sg2i

  function prep_ice_get_mapper_Fg2i()
    type(seq_map), pointer :: prep_ice_get_mapper_Fg2i
    prep_ice_get_mapper_Fg2i => mapper_Fg2i
  end function prep_ice_get_mapper_Fg2i

end module prep_ice_mod
