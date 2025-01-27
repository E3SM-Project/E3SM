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
  use seq_flds_mod, only: seq_flds_a2x_ext_fields, seq_flds_a2x_fields, seq_flds_x2a_fields ! 
  use seq_flds_mod, only: seq_flds_o2x_fields ! needed for MOAB init of ocean fields o2x to be able to transfer to coupler
  use seq_flds_mod, only: seq_flds_x2o_fields ! needed for MOAB init of ocean fields x2o to be able to transfer from coupler
  use seq_flds_mod, only: seq_flds_i2x_fields, seq_flds_x2i_fields ! needed for MOAB init of ice fields x2o on coupler side, to save them
  use seq_flds_mod, only: seq_flds_l2x_fields, seq_flds_x2l_fields ! 
  use seq_flds_mod, only: seq_flds_r2x_fields, seq_flds_x2r_fields !
  use seq_comm_mct, only: cplid, logunit
  use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs, seq_comm_iamin
  use seq_diag_mct

  use seq_comm_mct, only : mhid, mpoid, mbaxid, mboxid, mbofxid ! iMOAB app ids, for atm, ocean, ax mesh, ox mesh
  use seq_comm_mct, only : mhpgid         !    iMOAB app id for atm pgx grid, on atm pes
  use seq_comm_mct, only : atm_pg_active  ! flag if PG mesh instanced
  use seq_comm_mct, only : mlnid , mblxid !    iMOAB app id for land , on land pes and coupler pes
  use seq_comm_mct, only : mphaid !            iMOAB app id for phys atm; comp atm is 5, phys 5+200
  use seq_comm_mct, only : MPSIID, mbixid  !  sea-ice on comp pes and on coupler pes
  use seq_comm_mct, only : mrofid, mbrxid  ! iMOAB id of moab rof app on comp pes and on coupler too
  use shr_mpi_mod,  only: shr_mpi_max
  ! use dimensions_mod, only : np     ! for atmosphere
  use iso_c_binding

  implicit none
  private  ! except
#include <mpif.h>
#include "moab/MOABConfig.h"
  save

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: seq_map_init_exchange   ! union of cpl/component pes
  public :: seq_map_map_exchange    ! union of cpl/component pes
  public :: seq_mctext_gsmapInit
  public :: seq_mctext_avInit
  public :: seq_mctext_gGridInit
  public :: seq_mctext_avExtend
  public :: cplcomp_moab_Init       ! called to migrate MOAB mesh from
                                    !   component pes to coupler pes
  public :: component_exch_moab
  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  ! Shared routines for extension and computation of gsmaps, avs, and ggrids
  private :: seq_mctext_gsmapIdentical
  private :: seq_mctext_gsmapExtend
  private :: seq_mctext_gsmapCreate
  private :: seq_mctext_avCreate

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  integer,public :: seq_mctext_decomp

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(*),parameter :: subName = '(seq_mctext_mct)'
  real(r8),parameter :: c1 = 1.0_r8

  !=======================================================================
contains
  !=======================================================================

  subroutine seq_map_init_exchange( comp, mapper, flow, string)

    implicit none
    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(component_type), intent(inout)      :: comp
    type(seq_map)   , intent(inout), pointer :: mapper
    character(len=3), intent(in)             :: flow
    character(len=*), intent(in),optional    :: string
    !
    ! Local Variables
    !
    integer(IN)                :: ID_s
    integer(IN)                :: ID_d
    integer(IN)                :: ID_join
    integer(IN)                :: mapid
    integer(IN)                :: mpicom_s, mpicom_d, mpicom_join
    type(mct_gsmap) , pointer  :: gsmap_s
    type(mct_gsmap) , pointer  :: gsmap_d
    type(mct_gsmap)            :: gsmap_s_join
    type(mct_gsmap)            :: gsmap_d_join
    character(len=*),parameter :: subname = "(seq_map_init_rearrsplit) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    id_join = comp%cplcompid
    call seq_comm_getinfo(ID_join, mpicom=mpicom_join)

    if (flow == 'c2x') then
       gsmap_s => component_get_gsmap_cc(comp)
       gsmap_d => component_get_gsmap_cx(comp)
    end if
    if (flow == 'x2c') then
       gsmap_s => component_get_gsmap_cx(comp)
       gsmap_d => component_get_gsmap_cc(comp)
    end if

    if (mct_gsmap_Identical(gsmap_s,gsmap_d)) then

       call seq_map_mapmatch(mapid, gsmap_s=gsmap_s, gsmap_d=gsmap_d, strategy="copy")

       if (mapid > 0) then
          call seq_map_mappoint(mapid, mapper)
       else
          call seq_map_mapinit(mapper, mpicom_join)
          mapper%copy_only = .true.
          mapper%strategy = "copy"
          if (flow == 'c2x') then
             mapper%gsmap_s => component_get_gsmap_cc(comp)
             mapper%gsmap_d => component_get_gsmap_cx(comp)
          end if
          if (flow == 'x2c') then
             mapper%gsmap_s => component_get_gsmap_cx(comp)
             mapper%gsmap_d => component_get_gsmap_cc(comp)
          end if
       endif

       if (seq_comm_iamroot(ID_join)) then
          write(logunit,'(2A,L2)') subname,' gsmaps ARE IDENTICAL, copyoption = ',mapper%copy_only
       endif

    else

       if (seq_comm_iamroot(ID_join)) write(logunit,'(2A)') subname,' gsmaps are not identical'

       if (flow == 'c2x') then
          id_s = comp%compid
          id_d = cplid
       end if
       if (flow == 'x2c') then
          id_s = cplid
          id_d = comp%compid
       end if
       call seq_comm_getinfo(ID_s   , mpicom=mpicom_s)
       call seq_comm_getinfo(ID_d   , mpicom=mpicom_d)
       call seq_comm_getinfo(ID_join, mpicom=mpicom_join)

       ! --- Extend gsmaps to join group of pes

       call seq_mctext_gsmapExtend(gsmap_s, mpicom_s, gsmap_s_join, mpicom_join, ID_join)
       call seq_mctext_gsmapExtend(gsmap_d, mpicom_d, gsmap_d_join, mpicom_join, ID_join)

       ! --- Initialize rearranger based on join gsmaps
       ! --- test for the gsmaps instead of the gsmap joins because the gsmap joins are temporary

       ! -------------------------------
       ! tcx  tcraig mapmatch is a problem here because we're comparing gsmaps that may not be defined
       !      on some pes.  first issue is whether gsmap_identical in underlying routine will abort.
       !      second issue is whether different pes return different values.  use mapidmin, mapidmax to
       !      confirm all mapids returned are the same.  if not, then just set mapid to -1 and compute
       !      a new rearranger.
       ! tcx  not clear this works all the time, so just do not do map matching here for time being
       !      Sept 2013.
       ! -------------------------------
       !       mapid = -1
       !       call seq_map_mapmatch(mapid,gsmap_s=gsmap_s,gsmap_d=gsmap_d,strategy="rearrange")
       !       call shr_mpi_min(mapid,mapidmin,mpicom_join,subname//' min')
       !       call shr_mpi_max(mapid,mapidmax,mpicom_join,subname//' max')
       !       if (mapidmin /= mapidmax) mapid = -1
       ! -------------------------------

       ! --- Initialize rearranger
       ! --- the gsmap joins are temporary so store the regular gsmaps in the mapper
       call seq_map_mapinit(mapper, mpicom_join)
       mapper%rearrange_only = .true.
       mapper%strategy = "rearrange"
       if (flow == 'c2x') then
          mapper%gsmap_s => component_get_gsmap_cc(comp)
          mapper%gsmap_d => component_get_gsmap_cx(comp)
       end if
       if (flow == 'x2c') then
          mapper%gsmap_s => component_get_gsmap_cx(comp)
          mapper%gsmap_d => component_get_gsmap_cc(comp)
       end if
       call seq_map_gsmapcheck(gsmap_s_join, gsmap_d_join)
       call mct_rearr_init(gsmap_s_join, gsmap_d_join, mpicom_join, mapper%rearr)

       ! --- Clean up temporary gsmaps

       call mct_gsMap_clean(gsmap_s_join)
       call mct_gsMap_clean(gsmap_d_join)

    endif

    if (seq_comm_iamroot(CPLID)) then
       write(logunit,'(2A,I6,4A)') subname,' mapper counter, strategy, mapfile = ', &
            mapper%counter,' ',trim(mapper%strategy),' ',trim(mapper%mapfile)
       call shr_sys_flush(logunit)
    endif

  end subroutine seq_map_init_exchange

  !===============================================================================

  subroutine seq_map_map_exchange( comp, flow, dom_flag, dom_tmp, string, msgtag )

    !-----------------------------------------------------
    !
    ! Arguments
    !
    type(component_type) , intent(inout)               :: comp
    character(len=3)     , intent(in)                  :: flow
    logical              , intent(in),optional         :: dom_flag
    type(mct_gGrid)      , intent(in),optional, target :: dom_tmp
    character(len=*)     , intent(in),optional         :: string
    integer(IN)          , intent(in),optional         :: msgtag
    !
    ! Local Variables
    !
    type(seq_map)  , pointer :: mapper
    type(mct_aVect), pointer :: av_s
    type(mct_aVect), pointer :: av_d
    type(mct_gGrid), pointer :: dom_s
    type(mct_gGrid), pointer :: dom_d
    integer(IN),save         :: ltag    ! message tag for rearrange
    character(len=*),parameter :: subname = "(seq_map_map) "
    !-----------------------------------------------------

    if (seq_comm_iamroot(CPLID) .and. present(string)) then
       write(logunit,'(A)') subname//' called for '//trim(string)
    endif

    if (flow == 'c2x') then
       if (present(dom_flag)) then
          dom_s   => component_get_dom_cc(comp)
          dom_d   => component_get_dom_cx(comp)
          ! Overwrite dom_d pointer if dom_tmp is present
          ! Needed for backwards compatibility with domain checker in component_init_cx
          if (present(dom_tmp)) then
             dom_d => dom_tmp
          end if
       else
          av_s   => component_get_c2x_cc(comp)
          av_d   => component_get_c2x_cx(comp)
       end if
       mapper => component_get_mapper_Cc2x(comp)
    end if
    if (flow == 'x2c') then
       if (present(dom_flag)) then
          dom_s  => component_get_dom_cx(comp)
          dom_d  => component_get_dom_cc(comp)
       else
          av_s   => component_get_x2c_cx(comp)
          av_d   => component_get_x2c_cc(comp)
       end if
       mapper => component_get_mapper_Cx2c(comp)
    end if

    if (present(msgtag)) then
       ltag = msgtag
    else
       ltag = 2000
    endif

    if (mapper%copy_only) then
       !-------------------------------------------
       ! COPY data
       !-------------------------------------------
       if (present(dom_flag)) then
          call mct_aVect_copy(aVin=dom_s%data, aVout=dom_d%data, vector=mct_usevector)
       else
          call mct_aVect_copy(aVin=av_s, aVout=av_d, vector=mct_usevector)
       end if

    else if (mapper%rearrange_only) then
       !-------------------------------------------
       ! REARRANGE data
       !-------------------------------------------
       if (present(dom_flag)) then
          call mct_rearr_rearrange(dom_s%data, dom_d%data, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall)
       else
          call mct_rearr_rearrange(av_s, av_d, mapper%rearr, tag=ltag, VECTOR=mct_usevector, &
               ALLTOALL=mct_usealltoall)
       end if
    end if

  end subroutine seq_map_map_exchange

  !=======================================================================

  subroutine seq_mctext_gsmapInit(comp)

    ! This routine initializes a gsmap based on another gsmap potentially
    ! on other pes.  It addresses non-overlap of pes.

    !-----------------------------------------------------
    !
    ! Arguments
    !
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
    type(mct_gsMap), pointer :: gsmap_old
    type(mct_gsMap), pointer :: gsmap_new
    type(mct_gsMap)          :: gsmap_old_join   ! gsmap_old on joined id, temporary
    character(len=*),parameter :: subname = "(seq_mctext_gsmapInit) "
    !-----------------------------------------------------

    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)

    id_new  = CPLID
    id_old  = comp%compid
    id_join = comp%cplcompid

    mpicom_new  = mpicom_CPLID
    mpicom_old  = comp%mpicom_compid
    mpicom_join = comp%mpicom_cplcompid

    gsmap_new => component_get_gsmap_cx(comp)
    gsmap_old => component_get_gsmap_cc(comp)

    call seq_comm_getinfo(ID_old ,mpicom=mpicom_old)
    call seq_comm_getinfo(ID_new ,mpicom=mpicom_new)
    call seq_comm_getinfo(ID_join,mpicom=mpicom_join)

    ! --- Set gsmaps
    ! ---   Extend the old one to now span all pes on ID_join
    ! ---   Create a new gsmap on pes associated with ID_new using info from the old one

    call seq_mctext_gsmapExtend(gsmap_old     , mpicom_old  , gsmap_old_join, mpicom_join, ID_join)
    call seq_mctext_gsmapCreate(gsmap_old_join, mpicom_join , gsmap_new     , mpicom_new , ID_new  )

    call mct_gsMap_clean(gsmap_old_join)

  end subroutine seq_mctext_gsmapInit

  !=======================================================================

  subroutine seq_mctext_avInit( comp, flow )

    !-----------------------------------------------------
    ! This routine initializes Avs that may need to be extended
    !
    ! Arguments
    !
    type(component_type), intent(inout) :: comp
    character(len=3)    , intent(in)    :: flow
    !
    ! Local Variables
    !
    integer                  :: lsize
    integer                  :: mpicom_cplid
    integer                  :: mpicom_new
    integer                  :: ID_old
    integer                  :: ID_new
    integer                  :: ID_join
    type(mct_aVect), pointer :: AV1_old
    type(mct_aVect), pointer :: AV1_new
    type(mct_gsmap), pointer :: gsmap_new
    character(len=*),parameter :: subname = "(seq_mctext_avInit) "
    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_comm_getinfo(CPLID ,mpicom=mpicom_CPLID)

    id_new  = cplid
    id_old  = comp%compid
    id_join = comp%cplcompid

    mpicom_new  = mpicom_cplid

    gsmap_new => component_get_gsmap_cx(comp)

    if (flow == 'c2x') then
       av1_old => component_get_c2x_cc(comp)
       av1_new => component_get_c2x_cx(comp)
    end if
    if (flow == 'x2c') then
       av1_old => component_get_x2c_cc(comp)
       av1_new => component_get_x2c_cx(comp)
    end if

    ! --- Extend old avs and initialize new avs for use in the future

    lsize = 0
    if (seq_comm_iamin(ID_new)) then
       lsize = mct_gsMap_lsize(gsMap_new, mpicom_new)
    endif
    call seq_mctext_avExtend(AV1_old, ID_old, ID_join)
    call seq_mctext_avCreate(AV1_old, ID_old, AV1_new, ID_join, lsize)

  end subroutine seq_mctext_avInit

  !=======================================================================

  subroutine seq_mctext_gGridInit(comp, ggrid_new)

    !-----------------------------------------------------
    ! This routine initializes gGrids that may need to be extended
    !
    ! Arguments
    !
    type(component_type), intent(inout) :: comp
    type(mct_gGrid), optional, target, intent(inout) :: ggrid_new
    !
    ! Local Variables
    !
    integer                  :: mpicom_cplid
    integer                  :: lsize
    integer                  :: mpicom_new
    integer                  :: ID_old
    integer                  :: ID_new
    integer                  :: ID_join
    type(mct_gGrid), pointer :: GG1_old
    type(mct_gGrid), pointer :: GG1_new
    type(mct_gsmap), pointer :: gsmap_new
    character(len=*),parameter :: subname = "(seq_mctext_gGridInit) "
    !-----------------------------------------------------

    ! --- Setup data for use and make sure the old ID is ok

    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)

    id_new  = cplid
    id_old  = comp%compid
    id_join = comp%cplcompid

    mpicom_new = mpicom_cplid

    gsmap_new => component_get_gsmap_cx(comp)

    gg1_old => component_get_dom_cc(comp)
    gg1_new => component_get_dom_cx(comp)

    ! --- Extend old ggrids and initialize new ggrids for use in the future

    lsize = 0
    if (seq_comm_iamin(ID_new)) then
       lsize = mct_gsMap_lsize(gsMap_new,mpicom_new)
    endif
    call seq_mctext_avExtend(GG1_old%data, ID_old, ID_join)

    if (present(ggrid_new)) then
       call mct_gGrid_init(GGrid=ggrid_new, CoordChars=seq_flds_dom_coord, OtherChars=seq_flds_dom_other, lsize=lsize )
       call mct_avect_zero(ggrid_new%data)
    else
       call mct_gGrid_init(GGrid=GG1_new, CoordChars=seq_flds_dom_coord, OtherChars=seq_flds_dom_other, lsize=lsize )
       call mct_avect_zero(GG1_new%data)
    end if

  end subroutine seq_mctext_gGridInit

  !=======================================================================

  subroutine seq_mctext_gsmapExtend(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !----------------------------------------------------------------
    ! Extend/Convert a gsmap from one mpicom to another mpicom that contains
    ! at least all the pes that gsmap uses, but with different ranks
    !----------------------------------------------------------------

    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_mctext_gsmapExtend) "
    integer :: n
    integer :: ngseg
    integer :: gsize
    integer :: msizei,msizeo
    integer :: mrank,mranko,mrankog   ! sets pe rank of root mpicomi pe in mpicomo
    integer :: mpigrpi,mpigrpo
    integer :: ierr
    integer, pointer :: pei(:),peo(:)
    integer, pointer :: start(:),length(:),peloc(:)

    mranko = -1

    ! --- create the new gsmap on the mpicomi root only

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_rank i')
       if (mrank == 0) then
          call mpi_comm_group(mpicomi,mpigrpi,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group i')
          call mpi_comm_group(mpicomo,mpigrpo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_group o')
          call mpi_comm_size(mpicomi,msizei,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size i')
          call mpi_comm_size(mpicomo,msizeo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_comm_size o')

          ! --- setup the translation of pe numbers from the old gsmap(mpicom)
          ! --- to the new one, pei -> peo

          allocate(pei(0:msizei-1),peo(0:msizei-1))
          do n = 0,msizei-1
             pei(n) = n
          enddo

          peo = -1
          call mpi_group_translate_ranks(mpigrpi,msizei,pei,mpigrpo,peo,ierr)
          call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_group_translate_ranks')

          do n = 0,msizei-1
             if (peo(n) < 0 .or. peo(n) > msizeo-1) then
                write(logunit,*) subname,' peo out of bounds ',peo(n),msizeo
                call shr_sys_abort()
             endif
          enddo

          mranko = peo(0)

          ! --- compute the new gsmap which has the same start and length values
          ! --- but peloc is now the mapping of pei to peo

          ngseg = gsmapi%ngseg
          gsize = gsmapi%gsize
          allocate(start(ngseg),length(ngseg),peloc(ngseg))
          do n = 1,ngseg
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
             peloc(n)  = peo(gsmapi%pe_loc(n))
          enddo

          ! --- initialize the gsmap on the root pe

          call mct_gsmap_init(gsmapo,compido,ngseg,gsize,start,length,peloc)

          deallocate(pei,peo,start,length,peloc)
       endif
    endif

    ! --- broadcast via allreduce the mpicomi root pe in mpicomo space
    ! --- mranko is -1 except on the root pe where is it peo of that pe

    call mpi_allreduce(mranko,mrankog,1,MPI_INTEGER,MPI_MAX,mpicomo,ierr)
    call shr_mpi_chkerr(ierr,subname//' gsm_cop mpi_allreduce max')

    ! --- broadcast the gsmap to all pes in mpicomo from mrankog

    call mct_gsmap_bcast(gsmapo, mrankog, mpicomo)

    ! tcx summarize decomp info
#if (1 == 0)
    write(logunit,*) trim(subname),'tcxa ',mpicomi,mpicomo
    call shr_sys_flush(logunit)
    call mpi_barrier(mpicomo,ierr)

    if (mpicomi /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mrank,ierr)
       write(logunit,*) 'tcxbi ',mrank
       if (mrank == 0) then
          write(logunit,*) 'tcxci ',gsmapi%ngseg,size(gsmapi%start),gsmapi%gsize,gsmapi%comp_id
          do n = 1,gsmapi%ngseg
             write(logunit,*) 'tcx gsmti ',n,gsmapi%start(n),gsmapi%length(n),gsmapi%pe_loc(n)
          enddo
          call shr_sys_flush(logunit)
       endif
    endif

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomo,mrank,ierr)
       write(logunit,*) 'tcxbo ',mrank
       if (mrank == 0) then
          write(logunit,*) 'tcxco ',gsmapo%ngseg,size(gsmapo%start),gsmapo%gsize,gsmapo%comp_id
          do n = 1,gsmapo%ngseg
             write(logunit,*) 'tcx gsmto ',n,gsmapo%start(n),gsmapo%length(n),gsmapo%pe_loc(n)
          enddo
          call shr_sys_flush(logunit)
       endif
    endif

    call shr_sys_flush(logunit)
    call mpi_barrier(mpicomo,ierr)
#endif


  end subroutine seq_mctext_gsmapExtend

  !=======================================================================

  subroutine seq_mctext_gsmapCreate(gsmapi, mpicomi, gsmapo, mpicomo, compido)

    !---------------------------------------------------------------------
    ! creates a new gsmap on a subset of pes, requires setting a new decomp
    !---------------------------------------------------------------------

    implicit none
    type(mct_gsMap), intent(IN) :: gsmapi
    integer        , intent(IN) :: mpicomi
    type(mct_gsMap), intent(OUT):: gsmapo
    integer        , intent(IN) :: mpicomo
    integer        , intent(IN) :: compido

    character(len=*),parameter :: subname = "(seq_mctext_gsmapCreate) "
    integer :: n,m,k
    integer :: ktot            ! number of active cells in gsmap
    integer :: apesi, apeso    ! number of active pes in gsmap
    integer ::        lsizeo   ! local size for lindex
    integer :: ngsegi,ngsego   ! ngseg of mpicomi, mpicomo
    integer :: gsizei,gsizeo   ! gsize of mpicomi, mpicomo
    integer :: msizei,msizeo   ! size of mpicomi, mpicomo
    integer :: mranki,mranko   ! rank in mpicomi, mpicomo
    integer :: ierr
    integer :: decomp_type
    integer, pointer :: start(:),length(:),peloc(:),perm(:),gindex(:),lindex(:)
    real(r8):: rpeloc
    logical :: gsmap_bfbflag = .false. ! normally this should be set to false

    ! --- create a new gsmap on new pes based on the old gsmap
    ! --- gsmapi must be known on all mpicomo pes, compute the same
    ! --- thing on all pes in parallel

    if (mpicomo /= MPI_COMM_NULL) then
       call mpi_comm_rank(mpicomi,mranki,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank i')
       call mpi_comm_size(mpicomi,msizei,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size i')
       call mpi_comm_rank(mpicomo,mranko,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank o')
       call mpi_comm_size(mpicomo,msizeo,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size o')

       ngsegi = gsmapi%ngseg
       gsizei = gsmapi%gsize
       gsizeo = gsizei
       call mct_gsMap_activepes(gsmapi,apesi)

       decomp_type = 0

       if (seq_mctext_decomp == 0) then
          if (msizeo == apesi) then      ! preserve segments and decomp
             ! For testing - set decomp_type to 1 - to have gsmapi and gsmapo identical
             if (gsmap_bfbflag) then
                decomp_type = 1     ! better in cpl to have all decomps "same-ish"
             else
                decomp_type = 2
             end if
          elseif (ngsegi >= msizeo) then ! preserve segments, new decomp
             decomp_type = 2
          else                           ! new segments
             decomp_type = 3
          endif
       else
          decomp_type = seq_mctext_decomp
       endif

       !tcx       decomp_type = 3 ! over ride setting above for testing
       !       if (mranko == 0) write(logunit,'(2A,4I)') trim(subname),' decomp_type =',decomp_type,ngsegi,msizeo,apesi

       select case (decomp_type)

       case(1)   ! --- preserve segments and decomp ---------------------

          ! -- copy the gsmap and translate the pes
          call mct_gsMap_copy(gsmapi,gsmapo)
          ngsego = ngsegi
          do n = 1,ngsego
             gsmapo%pe_loc(n) = mod(gsmapo%pe_loc(n),msizeo)    ! translate pes 1:1 from old to new
          enddo

       case(2)   ! --- preserve segments, new decomp --------------------

          ! --- preserve segments, sort the start and length, assign a new pe list
          ngsego = ngsegi
          allocate(start(ngsego),length(ngsego),peloc(ngsego),perm(ngsego))
          do n = 1,ngsego
             start(n)  = gsmapi%start(n)
             length(n) = gsmapi%length(n)
          enddo
          ! --- sort gsmap to minimize permute cost in mct
          call mct_indexset(perm)
          call mct_indexsort(ngsego,perm,start)
          call mct_permute(start,perm,ngsego)
          call mct_permute(length,perm,ngsego)
          ! --- give each pe "equal" number of segments, use reals to avoid integer overflow
          do n = 1,ngsego
             rpeloc = (((msizeo*c1)*((n-1)*c1))/(ngsego*c1))      ! give each pe "equal" number of segments, use reals to avoid integer overflow
             peloc(n) = int(rpeloc)
          enddo
          call mct_gsmap_init(gsmapo,ngsego,start,length,peloc,0,mpicomo,compido,gsizeo)
          deallocate(start,length,peloc,perm)

       case(3)   ! --- new segments, new decomp -------------------------

          ! --- new segments, compute gindex, then parse the gridcells out evenly

          k = 0
          do n = 1,ngsegi
             do m = 1,gsmapi%length(n)
                k = k + 1
                if (k > gsizei) then
                   write(logunit,*) trim(subname),' ERROR in gindex ',k,gsizei
                   call shr_sys_abort()
                endif
             enddo
          enddo
          ktot = k

          allocate(gindex(ktot),perm(ktot))

          k = 0
          do n = 1,ngsegi
             do m = 1,gsmapi%length(n)
                k = k + 1
                gindex(k) = gsmapi%start(n) + m - 1
             enddo
          enddo
          call mct_indexset(perm)
          call mct_indexsort(ktot,perm,gindex)
          call mct_permute(gindex,perm,ktot)

          k = 0
          do m = 0,msizeo-1
             lsizeo = ktot/msizeo
             if (m < (ktot - lsizeo*msizeo)) lsizeo = lsizeo + 1
             if (mranko == m) then
                allocate(lindex(lsizeo))
                if (k+lsizeo > ktot) then
                   write(logunit,*) trim(subname),' ERROR: decomp out of bounds ',mranko,k,lsizeo,ktot
                   call shr_sys_abort()
                endif
                lindex(1:lsizeo) = gindex(k+1:k+lsizeo)
                !                write(logunit,*) trim(subname),' decomp is ',mranko,lsizeo,k+1,k+lsizeo
             endif
             k = k + lsizeo
          enddo
          if (k /= ktot) then
             write(logunit,*) trim(subname),' ERROR: decomp incomplete ',k,ktot
             call shr_sys_abort()
          endif

          call mct_gsmap_init(gsmapo,lindex,mpicomo,compido,size(lindex),gsizeo)
          deallocate(gindex,perm,lindex)

       case default   ! --- unknown ---
          write(logunit,*) trim(subname),' ERROR decomp_type unknown ',decomp_type
          call shr_sys_abort(trim(subname)//' ERROR decomp_type unknown')

       end select

       if (mranko == 0) then
          write(logunit,102) trim(subname),' created new gsmap decomp_type =',decomp_type
          write(logunit,102) trim(subname),'   ngseg/gsize        = ', &
               mct_gsmap_ngseg(gsmapo),mct_gsmap_gsize(gsmapo)
          call mct_gsmap_activepes(gsmapo,apeso)
          write(logunit,102) trim(subname),'   mpisize/active_pes = ', &
               msizeo,apeso
          write(logunit,102) trim(subname),'   avg seg per pe/ape = ', &
               mct_gsmap_ngseg(gsmapo)/msizeo,mct_gsmap_ngseg(gsmapo)/apeso
          write(logunit,102) trim(subname),'   nlseg/maxnlsegs    = ', &
               mct_gsmap_nlseg(gsmapo,0),mct_gsmap_maxnlseg(gsmapo)
102       format(2A,2I8)
       endif

       !       if (.not. mct_gsmap_increasing(gsmapo) ) then
       !          write(logunit,*) trim(subname),' ERROR: gsmapo not increasing'
       !          call shr_sys_abort()
       !       endif

    endif

  end subroutine seq_mctext_gsmapCreate

  !=======================================================================

  subroutine seq_mctext_avExtend(AVin,IDin,ID)

    !-----------------------------------------------------------------------
    ! Extend an AV to a larger set of pes or
    ! Initialize an AV on another set of pes
    !
    ! Arguments
    !
    type(mct_aVect), intent(INOUT):: AVin
    integer         ,intent(IN)   :: IDin ! ID associated with AVin
    integer        , intent(IN)   :: ID   ! ID to initialize over
    !
    ! Local variables
    !
    character(len=*),parameter :: subname = "(seq_mctext_avExtend) "
    integer :: mpicom
    integer :: rank
    integer :: lsizei, lsizen
    integer :: srank,srankg
    integer :: ierr
    character(len=CXX) :: iList,rList
    !-----------------------------------------------------------------------

    call seq_comm_getinfo(ID,mpicom=mpicom,iam=rank)

    ! --- lsizen is the size of the newly initialized AV, zero is valid
    ! --- lsizei is -1 on any peszero on any pes where AV is not yet initialized

    lsizei = -1
    if (seq_comm_iamin(IDin)) lsizei = mct_aVect_lsize(AVin)
    lsizen = 0

    ! --- find a pe that already has AVin allocated, use MPI_MAX to do so
    ! --- set the pe and broadcast it to all other pes using mpi_allreduce

    srank = -1
    srankg = -1
    if (lsizei > 0) srank = rank

    call mpi_allreduce(srank,srankg,1,MPI_INTEGER,MPI_MAX,mpicom,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_allreduce max')

    if (srankg < 0) then
       write(logunit,*) subname,' WARNING AVin empty '
       return
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

    ! --- now allocate the AV on any pes where the orig size is zero.  those
    ! --- should be pes that either have no data and may have been allocated
    ! --- before (no harm in doing it again) or have never been allocated

    if (lsizei <= 0) then
       if(len_trim(iList) > 0 .and. len_trim(rList) > 0) then
          call mct_aVect_init(AVin,iList=iList,rList=rList,lsize=lsizen)
       elseif (len_trim(iList) > 0 .and. len_trim(rList) == 0) then
          call mct_aVect_init(AVin,iList=iList,lsize=lsizen)
       elseif (len_trim(iList) == 0 .and. len_trim(rList) > 0) then
          call mct_aVect_init(AVin,rList=rList,lsize=lsizen)
       endif
    endif

  end subroutine seq_mctext_avExtend

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

  logical function seq_mctext_gsmapIdentical(gsmap1,gsmap2)

    implicit none
    type(mct_gsMap), intent(IN):: gsmap1
    type(mct_gsMap), intent(IN):: gsmap2

    ! Local variables

    character(len=*),parameter :: subname = "(seq_mctext_gsmapIdentical) "
    integer :: n
    logical :: identical

    !-----------------------

    identical = .true.

    ! --- continue compare ---
    if (identical) then
       if (mct_gsMap_gsize(gsmap1) /= mct_gsMap_gsize(gsmap2)) identical = .false.
       if (mct_gsMap_ngseg(gsmap1) /= mct_gsMap_ngseg(gsmap2)) identical = .false.
    endif

    ! --- continue compare ---
    if (identical) then
       do n = 1,mct_gsMap_ngseg(gsmap1)
          if (gsmap1%start(n)  /= gsmap2%start(n) ) identical = .false.
          if (gsmap1%length(n) /= gsmap2%length(n)) identical = .false.
          if (gsmap1%pe_loc(n) /= gsmap2%pe_loc(n)) identical = .false.
       enddo
    endif

    seq_mctext_gsmapIdentical = identical

  end function seq_mctext_gsmapIdentical

  !=======================================================================

   subroutine cplcomp_moab_Init(infodata,comp)

      ! This routine initializes an iMOAB app on the coupler pes,
      !  corresponding to the component pes. It uses send/receive 
      !  from iMOAB to replicate the mesh on coupler pes

      !-----------------------------------------------------
      !
      use iMOAB, only: iMOAB_RegisterApplication, iMOAB_ReceiveMesh, iMOAB_SendMesh, &
      iMOAB_WriteMesh, iMOAB_DefineTagStorage, iMOAB_GetMeshInfo, &
      iMOAB_FreeSenderBuffers, iMOAB_ComputeCommGraph, iMOAB_LoadMesh
      !
      use seq_infodata_mod
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
      character*200            :: appname, outfile, wopts, ropts
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
      integer                  nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
      ! type(mct_list)           :: temp_list
      ! integer                  :: nfields, arrsize
      ! real(R8), allocatable, target :: values(:)


   !-----------------------------------------------------

      call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID, iam=rank)

      id_new  = CPLID
      id_old  = comp%compid
      id_join = comp%cplcompid

      mpicom_new  = mpicom_CPLID
      mpicom_old  = comp%mpicom_compid
      mpicom_join = comp%mpicom_cplcompid

      partMethod = 0 ! trivial partitioning
      context_id = -1 ! original sends/receives, so the context is -1
                     ! needed only to free send buffers
#ifdef MOAB_HAVE_ZOLTAN
      partMethod = 2 ! it is better to use RCB for atmosphere and ocean too
#endif

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

!!!!!!!!!!!!!!!! ATMOSPHERE
      if ( comp%oneletterid == 'a' .and. maxMH /= -1) then
         call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
         call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

         ! find atm mesh/domain file if it exists; it would be for data atm model (atm_prognostic false)
         call seq_infodata_GetData(infodata,atm_mesh = atm_mesh)
         ! now, if on coupler pes, receive mesh; if on comp pes, send mesh

         if (mphaid >= 0) then  ! component atm procs
            ierr  = iMOAB_GetMeshInfo ( mphaid, nvert, nvise, nbl, nsurf, nvisBC )
            comp%mbApCCid = mphaid ! phys atm 
            comp%mbGridType = 0 ! point cloud
            comp%mblsize = nvert(1) ! point cloud
         endif

         if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (atmosphere)
         !  send mesh to coupler
            if ( trim(atm_mesh) == 'none' ) then ! full model
               if (atm_pg_active) then !  change : send the pg2 mesh, not coarse mesh, when atm pg active
                  ierr = iMOAB_SendMesh(mhpgid, mpicom_join, mpigrp_cplid, id_join, partMethod)
               else
                  ! still use the mhid, original coarse mesh
                  ierr = iMOAB_SendMesh(mhid, mpicom_join, mpigrp_cplid, id_join, partMethod)
               endif
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in sending mesh from atm comp '
                  call shr_sys_abort(subname//' ERROR in sending mesh from atm comp')
               endif
            endif
         endif ! atmosphere pes
         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            appname = "COUPLE_ATM"//C_NULL_CHAR
            ! migrated mesh gets another app id, moab atm to coupler (mbax)
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_new, id_join, mbaxid)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in registering ', appname
               call shr_sys_abort(subname//' ERROR registering '// appname)
            endif
            if ( trim(atm_mesh) == 'none' ) then ! full atm
               ierr = iMOAB_ReceiveMesh(mbaxid, mpicom_join, mpigrp_old, id_old)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in receiving mesh on atm coupler '
                  call shr_sys_abort(subname//' ERROR in receiving mesh on atm coupler ')
               endif
            else   ! data atm
              ! we need to read the atm mesh on coupler, from domain file 
               ierr = iMOAB_LoadMesh(mbaxid, trim(atm_mesh)//C_NULL_CHAR, &
                "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION;NO_CULLING", 0)
               if ( ierr /= 0 ) then
                  write(logunit,*) 'Failed to load atm domain mesh on coupler'
                  call shr_sys_abort(subname//' ERROR Failed to load atm domain mesh on coupler  ')
               endif
               if (seq_comm_iamroot(CPLID)) then
                  write(logunit,'(A)') subname//' load atm domain mesh from file '//trim(atm_mesh)
               endif
               ! right now, turn atm_pg_active to true
               atm_pg_active = .true. ! FIXME TODO 
               ! need to add global id tag to the app, it will be used in restart
               tagtype = 0  ! dense, integer
               numco = 1
               tagname='GLOBAL_ID'//C_NULL_CHAR
               ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in adding global id tag to atmx'
                  call shr_sys_abort(subname//' ERROR in adding global id tag to atmx ')
               endif
            endif
            
         endif  ! on coupler pes
         !  iMOAB_FreeSenderBuffers needs to be called after receiving the mesh

         if (mhid .ge. 0) then  ! we are on component atm pes
            if ( trim(atm_mesh) == 'none' ) then  ! full atmosphere
               context_id = id_join
               if (atm_pg_active) then! we send mesh from mhpgid app
                  ierr = iMOAB_FreeSenderBuffers(mhpgid, context_id)
               else
                  ierr = iMOAB_FreeSenderBuffers(mhid, context_id)
               endif
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in freeing send buffers '
                  call shr_sys_abort(subname//' ERROR in freeing send buffers')
               endif
            endif
         endif  ! component atm pes

         ! graph between atm phys, mphaid, and atm dyn on coupler, mbaxid
         ! phys atm group is mpigrp_old, coupler group is mpigrp_cplid
         typeA = 2 ! point cloud for mphaid
         typeB = 1 ! spectral elements
         if (atm_pg_active) then
            typeB = 3 ! in this case, we will have cells associated with DOFs as GLOBAL_ID tag
         endif
         ATM_PHYS_CID = 200 + id_old ! 200 + 5 for atm, see line  969   ATM_PHYS = 200 + ATMID ! in
                                    !  components/cam/src/cpl/atm_comp_mct.F90
                                    !  components/data_comps/datm/src/atm_comp_mct.F90 ! line 177 !! 

         ierr = iMOAB_ComputeCommGraph( mphaid, mbaxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
             typeA, typeB, ATM_PHYS_CID, id_join) ! ID_JOIN is now 6 

         ! we can receive those tags only on coupler pes, when mbaxid exists
         ! we have to check that before we can define the tag
         if (mbaxid .ge. 0 ) then   !  coupler pes
            tagtype = 1  ! dense, double

            if (atm_pg_active) then
              tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
              numco = 1 !  usually 1 value per cell
            else ! this is not supported now, but leave it here
              tagname = trim(seq_flds_a2x_ext_fields)//C_NULL_CHAR ! MOAB versions of a2x for spectral
              numco = 16 ! np*np !  usually 16 values per cell, GLL points; should be 4 x 4 = 16
            endif

            ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags on atm on coupler '
               call shr_sys_abort(subname//' ERROR in defining tags ')
            endif

            tagname = trim(seq_flds_x2a_fields)//C_NULL_CHAR ! TODO should be also x2a_ext for spectral case 
            ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags seq_flds_x2a_fields on atm on coupler '
               call shr_sys_abort(subname//' ERROR in defining tags ')
            endif

            !add the normalization tag
            tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mbaxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on atm on coupler '
               call shr_sys_abort(subname//' ERROR in defining tags ')
            endif
            ! also, frac, area,  masks has to come from atm mphaid, not from domain file reader
            ! this is hard to digest :(
            tagname = 'lat:lon:area:frac:mask'//C_NULL_CHAR
            call component_exch_moab(comp, mphaid, mbaxid, 0, tagname)

         endif ! coupler pes

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
      endif  ! comp%oneletterid == 'a'

!!!!!!!!!!!!!!!! OCEAN
      ! ocean
      if (comp%oneletterid == 'o'  .and. maxMPO /= -1) then
         call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
         call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

         ! find ocean domain file if it exists; it would be for data ocean model (ocn_prognostic false)
         call seq_infodata_GetData(infodata,ocn_domain=ocn_domain)

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
            if ( trim(ocn_domain) == 'none' ) then
               !  send mesh to coupler
               ierr = iMOAB_SendMesh(mpoid, mpicom_join, mpigrp_cplid, id_join, partMethod)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in sending ocean mesh to coupler '
                  call shr_sys_abort(subname//' ERROR in sending ocean mesh to coupler ')
               endif
               comp%mbGridType = 1 ! cells
               comp%mblsize = nvise(1) ! cells
            else
               comp%mbGridType = 0 ! vertices
               comp%mblsize = nvert(1) ! vertices
            endif
         endif
         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            appname = "COUPLE_MPASO"//C_NULL_CHAR
            ! migrated mesh gets another app id, moab ocean to coupler (mbox)
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_new, id_join, mboxid)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' cannot register ocn app on coupler '
               call shr_sys_abort(subname//' ERROR cannot register ocn app on coupler  ')
            endif
            if ( trim(ocn_domain) == 'none' ) then
               ierr = iMOAB_ReceiveMesh(mboxid, mpicom_join, mpigrp_old, id_old)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' cannot receive ocn mesh on coupler '
                  call shr_sys_abort(subname//' ERROR cannot receive ocn mesh on coupler  ')
               endif
            else
              ! we need to read the ocean mesh on coupler, from domain file 
               ierr = iMOAB_LoadMesh(mboxid, trim(ocn_domain)//C_NULL_CHAR, &
                "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;NO_CULLING;REPARTITION", 0)
               if ( ierr /= 0 ) then
                  write(logunit,*) 'Failed to load ocean domain mesh on coupler'
                  call shr_sys_abort(subname//' ERROR Failed to load ocean domain mesh on coupler  ')
               endif
               if (seq_comm_iamroot(CPLID)) then
                  write(logunit,'(A)') subname//' load ocn domain mesh from file '//trim(ocn_domain)
               endif
               ! need to add global id tag to the app, it will be used in restart
               tagtype = 0  ! dense, integer
               numco = 1
               tagname='GLOBAL_ID'//C_NULL_CHAR
               ierr = iMOAB_DefineTagStorage(mboxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in adding global id tag to ocnx'
                  call shr_sys_abort(subname//' ERROR in adding global id tag to ocnx ')
               endif
            endif  ! end of defining couplers copy of ocean mesh

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
         if (mpoid .ge. 0) then  ! we are on component ocn pes
            if ( trim(ocn_domain) == 'none' ) then
               context_id = id_join
               ierr = iMOAB_FreeSenderBuffers(mpoid, context_id)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in freeing buffers '
                  call shr_sys_abort(subname//' ERROR in freeing buffers ')
               endif
            endif
         endif
         ! in case of domain read, we need to compute the comm graph
         if ( trim(ocn_domain) /= 'none' ) then
            ! we are now on joint pes, compute comm graph between data ocn and coupler model ocn
            typeA = 2 ! point cloud on component PEs
            typeB = 3 ! full mesh on coupler pes, we just read it
            ierr = iMOAB_ComputeCommGraph( mpoid, mboxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
               typeA, typeB, id_old, id_join) 
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing comm graph for data ocn model '
               call shr_sys_abort(subname//' ERROR in computing comm graph for data ocn model ')
            endif
            ! also, frac, area,  masks has to come from ocean mpoid, not from domain file reader
            ! this is hard to digest :(
            tagname = 'area:frac:mask'//C_NULL_CHAR
            call component_exch_moab(comp, mpoid, mboxid, 0, tagname)
         endif 

         ! start copy
         ! do another ocean copy of the mesh on the coupler, just because So_fswpen field 
         ! would appear twice on original mboxid, once from xao states, once from o2x states
         id_join = id_join + 1000! kind of random 
         if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (ocean)
            if ( trim(ocn_domain) == 'none' ) then
               !  send mesh to coupler, the second time! a copy would be cheaper
               ierr = iMOAB_SendMesh(mpoid, mpicom_join, mpigrp_cplid, id_join, partMethod)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in sending ocean mesh to coupler the second time'
                  call shr_sys_abort(subname//' ERROR in sending ocean mesh to coupler the second time ')
               endif
            endif
         endif

         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            appname = "COUPLE_MPASOF"//C_NULL_CHAR
            ! migrated mesh gets another app id, moab ocean to coupler (mbox)
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_new, id_join, mbofxid)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' cant register second ocean mesh on coupler'
               call shr_sys_abort(subname//' ERROR cant register second ocean mesh on coupler')
            endif
            if ( trim(ocn_domain) == 'none' ) then
               ierr = iMOAB_ReceiveMesh(mbofxid, mpicom_join, mpigrp_old, id_old)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' cant receive second ocean mesh on coupler'
                  call shr_sys_abort(subname//' ERROR cant receive second ocean mesh on coupler')
               endif
            else
                ! we need to read the ocean mesh on coupler, from domain file 
               ierr = iMOAB_LoadMesh(mbofxid, trim(ocn_domain)//C_NULL_CHAR, &
                "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;NO_CULLING;REPARTITION", 0)
               if ( ierr /= 0 ) then
                  write(logunit,*) 'Failed to load second ocean domain mesh on coupler'
                  call shr_sys_abort(subname//' ERROR Failed to load second ocean domain mesh on coupler  ')
               endif
               if (seq_comm_iamroot(CPLID)) then
                  write(logunit,'(A)') subname//' load ocn domain mesh from file for second ocn instance '//trim(ocn_domain)
               endif
               ! need to add global id tag to the app, it will be used in restart
               tagtype = 0  ! dense, integer
               numco = 1
               tagname='GLOBAL_ID'//C_NULL_CHAR
               ierr = iMOAB_DefineTagStorage(mbofxid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in adding global id tag to ocnx'
                  call shr_sys_abort(subname//' ERROR in adding global id tag to ocnx ')
               endif
            endif

            tagtype = 1  ! dense, real
            numco = 1
            tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mbofxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on ocn on coupler '
               call shr_sys_abort(subname//' ERROR in defining tags ')
            endif

         endif

         if (mpoid .ge. 0) then  ! we are on component ocn pes again, release buffers 
             if ( trim(ocn_domain) == 'none' ) then
               context_id = id_join
               ierr = iMOAB_FreeSenderBuffers(mpoid, context_id)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in freeing buffers '
                  call shr_sys_abort(subname//' ERROR in freeing buffers ')
               endif
             endif
         endif
         ! end copy 
#ifdef MOABDEBUG
         outfile = 'recMeshOcnF.h5m'//C_NULL_CHAR
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
      !  write out the mesh file to disk
         ierr = iMOAB_WriteMesh(mbofxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
             write(logunit,*) subname,' error in writing ocean mesh coupler '
             call shr_sys_abort(subname//' ERROR in writing ocean mesh coupler ')
         endif
#endif
      endif  ! end ocean

!!!!!!!!!!!!!!!! LAND
   !   land
      if (comp%oneletterid == 'l'  .and. maxMLID /= -1) then
         call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
         call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

         ! use land full mesh 
         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            appname = "COUPLE_LAND"//C_NULL_CHAR
            ! migrated mesh gets another app id, moab land to coupler (mblx)
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_new, id_join, mblxid)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in registering coupler land '
               call shr_sys_abort(subname//' ERROR in registering coupler land')
            endif
            ! do not receive the mesh anymore, read it from file, then pair it with mlnid, component land PC mesh
            ! similar to rof mosart mesh  
            
            ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION'//C_NULL_CHAR
            call seq_infodata_GetData(infodata,lnd_domain=lnd_domain)
            outfile = trim(lnd_domain)//C_NULL_CHAR
            nghlay = 0 ! no ghost layers 
            if (seq_comm_iamroot(CPLID) ) then
               write(logunit, *) "load land domain file from file: ", trim(lnd_domain), &
                 " with options: ", trim(ropts)
            endif
            ierr = iMOAB_LoadMesh(mblxid, outfile, ropts, nghlay)
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in reading land coupler mesh from ', trim(lnd_domain)
               call shr_sys_abort(subname//' ERROR in reading land coupler mesh')
            endif
            if (seq_comm_iamroot(CPLID)) then
               write(logunit,'(A)') subname//' load lnd domain mesh from file '//trim(lnd_domain)
            endif
            ! need to add global id tag to the app, it will be used in restart
            tagtype = 0  ! dense, integer
            numco = 1
            tagname='GLOBAL_ID'//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in adding global id tag to lndx'
               call shr_sys_abort(subname//' ERROR in adding global id tag to lndx ')
            endif

   !  need to define tags on land too
            tagname = trim(seq_flds_l2x_fields)//C_NULL_CHAR 
            tagtype = 1  ! dense, double
            numco = 1 !  one value per cell
            ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags l2x on coupler land'
               call shr_sys_abort(subname//' ERROR in defining tags l2x on coupler ')
            endif
            ! need also to define seq_flds_x2l_fields on coupler instance, and on land comp instance
            tagname = trim(seq_flds_x2l_fields)//C_NULL_CHAR 
            ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags x2l on coupler land'
               call shr_sys_abort(subname//' ERROR in defining tags x2l on coupler land')
            endif        

            !add the normalization tag
            tagname = trim(seq_flds_dom_fields)//":norm8wt"//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mblxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in defining tags seq_flds_dom_fields on lnd on coupler '
               call shr_sys_abort(subname//' ERROR in defining tags ')
            endif

         endif ! end of coupler pes

         ! we are now on joint pes, compute comm graph between lnd and coupler model 
         typeA = 2 ! point cloud on component PEs, land
         typeB = 3 ! full mesh on coupler pes, we just read it
         if (mlnid >= 0) then
            ierr  = iMOAB_GetMeshInfo ( mlnid, nvert, nvise, nbl, nsurf, nvisBC )
            comp%mbApCCid = mlnid ! phys atm 
            comp%mbGridType = typeA - 2 ! 0 or 1, pc or cells 
            comp%mblsize = nvert(1) ! vertices
         endif
         ierr = iMOAB_ComputeCommGraph( mlnid, mblxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
             typeA, typeB, id_old, id_join) 
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in computing comm graph for lnd model '
            call shr_sys_abort(subname//' ERROR in computing comm graph for lnd model ')
         endif

         tagname = 'lat:lon:area:frac:mask'//C_NULL_CHAR
         call component_exch_moab(comp, mlnid, mblxid, 0, tagname)

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

      endif  ! End of land model

!!!!!!!!!!!!!!!! SEA ICE
      ! sea - ice
      if (comp%oneletterid == 'i'  .and. maxMSID /= -1) then
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
               comp%mbGridType = 1 ! 0 or 1, pc or cells 
               comp%mblsize = nvise(1) ! cells   
               !  send sea ice mesh to coupler
               ierr = iMOAB_SendMesh(MPSIID, mpicom_join, mpigrp_cplid, id_join, partMethod)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in sending sea ice mesh to coupler '
                  call shr_sys_abort(subname//' ERROR in sending sea ice mesh to coupler ')
               endif
            else
               comp%mbGridType = 0 ! 0 or 1, pc or cells
               comp%mblsize = nvert(1) ! vertices
            endif
         endif 
         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            appname = "COUPLE_MPASSI"//C_NULL_CHAR
            ! migrated mesh gets another app id, moab moab sea ice to coupler (mbix)
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_new, id_join, mbixid)
            if ( trim(ice_domain) == 'none' ) then ! regular ice model
               ierr = iMOAB_ReceiveMesh(mbixid, mpicom_join, mpigrp_old, id_old)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in receiving ice mesh in coupler '
                  call shr_sys_abort(subname//' ERROR in receiving sea ice mesh in coupler ')
               endif
            else
               ! we need to read the mesh ice (domain file) 
               ierr = iMOAB_LoadMesh(mbixid, trim(ice_domain)//C_NULL_CHAR, &
                "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;NO_CULLING;REPARTITION", 0)
               if ( ierr /= 0 ) then
                  write(logunit,*) 'Failed to load ice domain mesh on coupler'
                  call shr_sys_abort(subname//' ERROR Failed to load ice domain mesh on coupler  ')
               endif
               if (seq_comm_iamroot(CPLID)) then
                  write(logunit,'(A)') subname//' load ice domain mesh from file '//trim(ice_domain)
               endif
               ! need to add global id tag to the app, it will be used in restart
               tagtype = 0  ! dense, integer
               numco = 1
               tagname='GLOBAL_ID'//C_NULL_CHAR
               ierr = iMOAB_DefineTagStorage(mbixid, tagname, tagtype, numco,  tagindex )
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' error in adding global id tag to icex'
                  call shr_sys_abort(subname//' ERROR in adding global id tag to icex ')
               endif
            endif ! end data ice

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
         if (MPSIID .ge. 0) then  ! we are on component sea ice pes
            if ( trim(ice_domain) == 'none' ) then
               context_id = id_join
               ierr = iMOAB_FreeSenderBuffers(MPSIID, context_id)
               if (ierr .ne. 0) then
                 write(logunit,*) subname,' error in freeing buffers '
                 call shr_sys_abort(subname//' ERROR in freeing buffers ')
               endif
            endif
         endif
        ! in case of ice domain read, we need to compute the comm graph
        if ( trim(ice_domain) /= 'none' ) then
            ! we are now on joint pes, compute comm graph between data ice and coupler model ice
            typeA = 2 ! point cloud on component PEs
            typeB = 3 ! full mesh on coupler pes, we just read it
            ierr = iMOAB_ComputeCommGraph( MPSIID, mbixid, mpicom_join, mpigrp_old, mpigrp_cplid, &
               typeA, typeB, id_old, id_join) 
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in computing comm graph for data ice model '
               call shr_sys_abort(subname//' ERROR in computing comm graph for data ice model ')
            endif
            ! also, frac, area,  masks has to come from ice MPSIID , not from domain file reader
            ! this is hard to digest :(
            tagname = 'area:frac:mask'//C_NULL_CHAR
            call component_exch_moab(comp, MPSIID, mbixid, 0, tagname)
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

      endif

!!!!!!!!!!!!!!!! RIVER
     ! rof
      if (comp%oneletterid == 'r'  .and. maxMRID /= -1) then
         call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
         call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

         if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
            appname = "COUPLE_MROF"//C_NULL_CHAR
            ierr = iMOAB_RegisterApplication(trim(appname), mpicom_new, id_join, mbrxid)

            ! load mesh from scrip file passed from river model
            call seq_infodata_GetData(infodata,rof_mesh=rtm_mesh,rof_domain=rof_domain)
            if ( trim(rof_domain) == 'none' ) then
               outfile = trim(rtm_mesh)//C_NULL_CHAR
               ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN'//C_NULL_CHAR
            else
               outfile = trim(rof_domain)//C_NULL_CHAR
               ropts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION'//C_NULL_CHAR
            endif
            nghlay = 0 ! no ghost layers 
            ierr = iMOAB_LoadMesh(mbrxid, outfile, ropts, nghlay)
            if (seq_comm_iamroot(CPLID)) then
               write(logunit,'(A)') subname//' load rof from file '//trim(outfile)
            endif
            if ( ierr .ne. 0  ) then
               call shr_sys_abort( subname//' ERROR: cannot read rof mesh on coupler' )
            end if
            
             ! need to add global id tag to the app, it will be used in restart
            tagtype = 0  ! dense, integer
            numco = 1
            tagname='GLOBAL_ID'//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
            if (ierr .ne. 0) then
               write(logunit,*) subname,' error in adding global id tag to rof'
               call shr_sys_abort(subname//' ERROR in adding global id tag to rof ')
            endif
            
            tagtype = 1  ! dense, double
            numco = 1 !  one value per cell / entity
            tagname = trim(seq_flds_r2x_fields)//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
            if ( ierr == 1 ) then
               call shr_sys_abort( subname//' ERROR: cannot define tags for rof on coupler' )
            end if
            tagname = trim(seq_flds_x2r_fields)//C_NULL_CHAR
            ierr = iMOAB_DefineTagStorage(mbrxid, tagname, tagtype, numco,  tagindex )
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

         if (mrofid >= 0) then  ! component pes
            ierr  = iMOAB_GetMeshInfo ( mrofid, nvert, nvise, nbl, nsurf, nvisBC )
            comp%mbApCCid = mrofid ! 
            comp%mbGridType = 0 ! 0 or 1, pc or cells 
            comp%mblsize = nvert(1) ! vertices
         endif

         ! we are now on joint pes, compute comm graph between rof and coupler model 
         typeA = 2 ! point cloud on component PEs
         typeB = 3 ! full mesh on coupler pes, we just read it
         ierr = iMOAB_ComputeCommGraph( mrofid, mbrxid, mpicom_join, mpigrp_old, mpigrp_cplid, &
             typeA, typeB, id_old, id_join) 
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in computing comm graph for rof model '
            call shr_sys_abort(subname//' ERROR in computing comm graph for rof model ')
         endif

         tagname = 'area:lon:lat:frac:mask'//C_NULL_CHAR
         call component_exch_moab(comp, mrofid, mbrxid, 0, tagname)

         ! if (mrofid .ge. 0) then  ! we are on component rof  pes
         !    context_id = id_join
         !    ierr = iMOAB_FreeSenderBuffers(mrofid, context_id)
         !    if (ierr .ne. 0) then
         !    write(logunit,*) subname,' error in freeing buffers '
         !    call shr_sys_abort(subname//' ERROR in freeing buffers ')
         !    endif
         ! endif
#ifdef MOABDEBUG
         outfile = 'recMeshRof.h5m'//C_NULL_CHAR
         wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
  !      write out the mesh file to disk
         ierr = iMOAB_WriteMesh(mbrxid, trim(outfile), trim(wopts))
         if (ierr .ne. 0) then
            write(logunit,*) subname,' error in writing rof mesh on coupler '
            call shr_sys_abort(subname//' ERROR in writing rof mesh on coupler ')
         endif
#endif
      endif ! end for rof coupler set up 

   end subroutine cplcomp_moab_Init


  ! can exchange data between mesh in component and mesh on coupler.  Either way.
  ! used in first hop of 2-hop
  subroutine component_exch_moab(comp, mbAPPid1, mbAppid2, direction, fields )

   use iMOAB ,  only: iMOAB_SendElementTag, iMOAB_ReceiveElementTag, iMOAB_WriteMesh, iMOAB_FreeSenderBuffers
   use seq_comm_mct, only :  num_moab_exports ! for debugging
   use ISO_C_BINDING, only : C_NULL_CHAR
   use shr_kind_mod      , only :  CXX => shr_kind_CXX
   !---------------------------------------------------------------
    ! Description
    ! send tags (fields) from component to coupler or from coupler to component

    type(component_type)     , intent(in)           :: comp
    ! direction 0 is from component to coupler; 1 is from coupler to component
    integer,                   intent(in)           :: mbAPPid1, mbAppid2, direction
    character(CXX)           , intent(in)           :: fields

    character(*), parameter :: subname = '(component_exch_moab)'
    integer :: id_join, source_id, target_id, ierr
    integer :: mpicom_join
    character(CXX)              :: tagname
    character*100 outfile, wopts, lnum, dir

  ! how to get mpicomm for joint comp + coupler
    id_join = comp%cplcompid
    call seq_comm_getinfo(ID_join,mpicom=mpicom_join)
!
    tagName = trim(fields)//C_NULL_CHAR

    if (direction .eq. 0) then
       source_id = comp%compid
       target_id = comp%cplcompid
    else ! direction eq 1
       source_id = comp%cplcompid
       target_id = comp%compid
    endif
    ! for atm, add 200 to component side, because we will involve always the point cloud
    ! we are not supporting anymore the spectral case, at least for the time being
    ! we need to fix fv-cgll projection first
    if (comp%oneletterid == 'a' .and. direction .eq. 0 ) then
       source_id = source_id + 200
    endif
    if (comp%oneletterid == 'a' .and. direction .eq. 1 ) then
       target_id = target_id + 200
    endif
    if (mbAPPid1 .ge. 0) then !  send

       ! basically, use the initial partitioning
       ierr = iMOAB_SendElementTag(mbAPPid1, tagName, mpicom_join, target_id)
       if (ierr .ne. 0) then
          call shr_sys_abort(subname//' cannot send element tag: '//trim(tagName))
       endif

    endif
    if ( mbAPPid2 .ge. 0 ) then !  we are on receiving end
       ierr = iMOAB_ReceiveElementTag(mbAPPid2, tagName, mpicom_join, source_id)
       if (ierr .ne. 0) then
          call shr_sys_abort(subname//' cannot receive element tag: '//trim(tagName))
       endif
    endif

!     ! we can now free the sender buffers
    if (mbAPPid1 .ge. 0) then
       ierr = iMOAB_FreeSenderBuffers(mbAPPid1, target_id)
       if (ierr .ne. 0) then
          call shr_sys_abort(subname//' cannot free sender buffers')
       endif
    endif
    
#ifdef MOABDEBUG
    if (direction .eq. 0 ) then
       dir = 'c2x'
    else
       dir = 'x2c'
    endif
    if (seq_comm_iamroot(CPLID) ) then
       write(logunit,'(A)') subname//' '//comp%ntype//' called in direction '//trim(dir)//' for fields '//trim(tagname)
    endif
    if (mbAPPid2 .ge. 0 ) then !  we are on receiving pes, for sure
      ! number_proj = number_proj+1 ! count the number of projections
      write(lnum,"(I0.2)") num_moab_exports
      
      outfile = comp%ntype//'_'//trim(dir)//'_'//trim(lnum)//'.h5m'//C_NULL_CHAR
      wopts   = ';PARALLEL=WRITE_PART'//C_NULL_CHAR !
      ierr = iMOAB_WriteMesh(mbAPPid2, trim(outfile), trim(wopts))
      if (ierr .ne. 0) then
          call shr_sys_abort(subname//' cannot write file '// outfile)
       endif
    endif
#endif

  end subroutine component_exch_moab

end module cplcomp_exchange_mod
