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
  use seq_comm_mct, only: cplid, logunit
  use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs, seq_comm_iamin
  use seq_diag_mct

  use seq_comm_mct, only : mhid, mpoid, mbaxid, mboxid  ! iMOAB app ids, for atm, ocean, ax mesh, ox mesh
  use shr_mpi_mod,  only: shr_mpi_max

  implicit none
  private  ! except
#include <mpif.h>
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

    id_new  = cplid
    id_old  = comp%compid
    id_join = comp%cplcompid

    mpicom_new  = mpicom_cplid
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

  subroutine cplcomp_moab_Init(comp)

    ! This routine initializes an iMOAB app on the coupler pes,
    !  corresponding to the component pes. It uses send/receive 
    !  from iMOAB to replicate the mesh on coupler pes

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

    character(len=*),parameter :: subname = "(cplcomp_moab_Init) "

    integer                  :: mpigrp_cplid ! coupler pes
    integer                  :: mpigrp_old   !  component group pes
    integer, external        :: iMOAB_RegisterFortranApplication, iMOAB_ReceiveMesh, iMOAB_SendMesh
    integer, external        :: iMOAB_WriteMesh, iMOAB_DefineTagStorage
    integer                  :: ierr
    character*32             :: appname, outfile, wopts, tagnameProj
    integer                  :: maxMH, maxMPO ! max pids for moab apps
    integer                  :: tagtype, numco,  tagindex

    !-----------------------------------------------------

    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID)

    id_new  = cplid
    id_old  = comp%compid
    id_join = comp%cplcompid

    mpicom_new  = mpicom_cplid
    mpicom_old  = comp%mpicom_compid
    mpicom_join = comp%mpicom_cplcompid

    call seq_comm_getinfo(ID_old ,mpicom=mpicom_old)
    call seq_comm_getinfo(ID_new ,mpicom=mpicom_new)
    call seq_comm_getinfo(ID_join,mpicom=mpicom_join)

    call shr_mpi_max(mhid, maxMH, mpicom_join, all=.true.)
    call shr_mpi_max(mpoid, maxMPO, mpicom_join, all=.true.)
    if (seq_comm_iamroot(CPLID) ) then
       write(logunit, *) "MOAB coupling:  maxMH: ", maxMH, " maxMPO: ", maxMPO
    endif
    ! this works now for atmosphere;
    if ( comp%oneletterid == 'a' .and. maxMH /= -1) then
      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes
      ! now, if on coupler pes, receive mesh; if on comp pes, send mesh
      if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (atmosphere)
        !  send mesh to coupler
        ierr = iMOAB_SendMesh(mhid, mpicom_join, mpigrp_cplid, id_join);
      endif
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
        appname = "COUPLE_ATM"//CHAR(0)
        ! migrated mesh gets another app id, moab atm to coupler (mbax)
        ierr = iMOAB_RegisterFortranApplication(trim(appname), mpicom_new, id_join, mbaxid)
        ierr = iMOAB_ReceiveMesh(mbaxid, mpicom_join, mpigrp_old, id_old)
        ! debug test
        outfile = 'recMeshAtm.h5m'//CHAR(0)
        wopts   = ';PARALLEL=WRITE_PART'//CHAR(0)
!      write out the mesh file to disk
        ierr = iMOAB_WriteMesh(mbaxid, trim(outfile), trim(wopts))
      endif
    endif
    if (comp%oneletterid == 'o'  .and. maxMPO /= -1) then
      call seq_comm_getinfo(cplid ,mpigrp=mpigrp_cplid)  ! receiver group
      call seq_comm_getinfo(id_old,mpigrp=mpigrp_old)   !  component group pes

      if (MPI_COMM_NULL /= mpicom_old ) then ! it means we are on the component pes (atmosphere)
        !  send mesh to coupler
        ierr = iMOAB_SendMesh(mpoid, mpicom_join, mpigrp_cplid, id_join);
      endif
      if (MPI_COMM_NULL /= mpicom_new ) then !  we are on the coupler pes
        appname = "COUPLE_MPASO"//CHAR(0)
        ! migrated mesh gets another app id, moab ocean to coupler (mbox)
        ierr = iMOAB_RegisterFortranApplication(trim(appname), mpicom_new, id_join, mboxid)
        ierr = iMOAB_ReceiveMesh(mboxid, mpicom_join, mpigrp_old, id_old)
        ! debug test
        outfile = 'recMeshOcn.h5m'//CHAR(0)
        wopts   = ';PARALLEL=WRITE_PART'//CHAR(0) !

        ! define here the tag that will be projected from atmosphere
        tagnameProj = 'a2oTAG_proj'//CHAR(0)
        tagtype = 1  ! dense, double
        numco = 1 !  one value per cell
        ierr = iMOAB_DefineTagStorage(mboxid, tagnameProj, tagtype, numco,  tagindex )

!      write out the mesh file to disk
        ierr = iMOAB_WriteMesh(mboxid, trim(outfile), trim(wopts))
      endif
    endif



  end subroutine cplcomp_moab_Init

end module cplcomp_exchange_mod
