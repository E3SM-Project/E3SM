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

end module component_type_mod
