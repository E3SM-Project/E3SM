module prep_aoflux_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8 
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_xao, num_inst_frc, num_inst_ocn
  use seq_comm_mct,     only: CPLID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs 
  use seq_infodata_mod, only: seq_infodata_getdata, seq_infodata_type   
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: atm, ocn

  implicit none
  private ! except
  save

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_aoflux_init

  public :: prep_aoflux_calc_xao_ox
  public :: prep_aoflux_calc_xao_ax

  public :: prep_aoflux_get_xao_ox
  public :: prep_aoflux_get_xao_ax
  public :: prep_aoflux_get_o2x_ax

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! attribute vectors
  type(mct_aVect), pointer :: xao_ox(:)   ! Atm-ocn fluxes, ocn grid, cpl pes 
  type(mct_aVect), pointer :: xao_ax(:)   ! Atm-ocn fluxes, atm grid, cpl pes 
  type(mct_aVect), pointer :: o2x_ax(:)   ! only needed if aoflux_grid is 'atm'

  ! seq_comm_getData variables
  logical :: iamroot_CPLID                ! .true. => CPLID masterproc
  integer :: mpicom_CPLID                 ! MPI cpl communicator

  ! seq_infodata_getData variables
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_aoflux_init (infodata, fractions_ox, fractions_ax)

    !---------------------------------------------------------------
    ! Description
    ! Initialize atm/ocn flux component and compute ocean albedos
    ! module variables
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    type(mct_aVect)          , intent(in)    :: fractions_ox(:) 
    type(mct_aVect)          , intent(in)    :: fractions_ax(:) 
    !
    ! Local Variables
    integer                     :: exi  , efi, eoi
    integer                     :: lsize_o
    integer                     :: lsize_a
    character(SHR_KIND_CS)      :: aoflux_grid ! grid for atm ocn flux calc
    type(mct_avect) , pointer   :: a2x_ax
    type(mct_avect) , pointer   :: o2x_ox
    character(*)    , parameter :: subname = '(prep_aoflux_init)'
    !---------------------------------------------------------------

    call seq_infodata_getdata(infodata,  &
         aoflux_grid=aoflux_grid)

    call seq_comm_getdata(CPLID, &
         mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

    a2x_ax => component_get_c2x_cx(atm(1)) 
    if (associated(a2x_ax)) then
       lsize_a = mct_aVect_lsize(a2x_ax)
    else
       lsize_a = 0
    end if

    o2x_ox => component_get_c2x_cx(ocn(1)) 
    if (associated(o2x_ox)) then 
       lsize_o = mct_aVect_lsize(o2x_ox)
    else
       lsize_o = 0
    end if

    allocate(xao_ax(num_inst_xao))
    do exi = 1,num_inst_xao
       call mct_aVect_init(xao_ax(exi), rList=seq_flds_xao_fields, lsize=lsize_a)
       call mct_aVect_zero(xao_ax(exi))
    end do
    allocate(xao_ox(num_inst_xao))
    do exi = 1,num_inst_xao
       call mct_aVect_init(xao_ox(exi), rList=seq_flds_xao_fields, lsize=lsize_o)
       call mct_aVect_zero(xao_ox(exi))
    enddo

    if (aoflux_grid == 'atm') then
       allocate(o2x_ax(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_ax(eoi), rList=seq_flds_o2x_fields, lsize=lsize_a)
          call mct_aVect_zero(o2x_ax(eoi))
       enddo
    end if

  end subroutine prep_aoflux_init

  !================================================================================================

  subroutine prep_aoflux_calc_xao_ax(fractions_ox, flds, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create xao_ox 
    !
    ! Uses
    use prep_atm_mod, only: prep_atm_get_mapper_So2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fo2a
    !
    ! Arguments
    type(mct_aVect) , intent(in)    :: fractions_ox(:)
    character(len=*), intent(in)    :: flds
    character(len=*), intent(in)    :: timer
    !
    ! Local Variables
    type(seq_map)   , pointer :: mapper_So2a
    type(seq_map)   , pointer :: mapper_Fo2a
    integer :: exi, efi
    character(*), parameter :: subname = '(prep_aoflux_calc_xao_ax)'
    character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (trim(flds) == 'albedos') then
       do exi = 1,num_inst_xao
          efi = mod((exi-1),num_inst_frc) + 1

          mapper_So2a => prep_atm_get_mapper_So2a()
          call seq_map_map(mapper_So2a, xao_ox(exi), xao_ax(exi), &
               fldlist=seq_flds_xao_albedo, norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
       enddo
    end if

    if (trim(flds) == 'states_and_fluxes') then
       do exi = 1,num_inst_xao
          efi = mod((exi-1),num_inst_frc) + 1

          mapper_So2a => prep_atm_get_mapper_So2a()
          call seq_map_map(mapper_So2a, xao_ox(exi), xao_ax(exi), &
               fldlist=seq_flds_xao_states, norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')

          mapper_Fo2a => prep_atm_get_mapper_Fo2a()
          call seq_map_map(mapper_Fo2a, xao_ox(exi), xao_ax(exi),&
               fldlist=seq_flds_xao_fluxes, norm=.true., &
               avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
       enddo
    end if
    call t_drvstopf  (trim(timer))

  end subroutine prep_aoflux_calc_xao_ax

  !================================================================================================

  subroutine prep_aoflux_calc_xao_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create xao_ox 
    !
    ! Uses
    use prep_ocn_mod, only: prep_ocn_get_mapper_Fa2o
    !
    ! Arguments
    character(len=*), intent(in)    :: timer
    !
    ! Local Variables
    type(seq_map), pointer :: mapper_Fa2o
    integer :: exi
    character(*), parameter :: subname = '(prep_aoflux_calc_xao_ax)'
    character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    ! this mapping has to be done with area overlap mapping for all fields 
    ! due to the masking of the xao_ax data and the fact that a2oS is bilinear

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do exi = 1,num_inst_xao 
       if (iamroot_CPLID .and. exi == 1) then 
          write(logunit,F00) 'Calling map_atm2ocn_mct for mapping xao_ax to xao_ox'
       end if

       mapper_Fa2o => prep_ocn_get_mapper_Fa2o()
       call seq_map_map(mapper_Fa2o, xao_ax(exi), xao_ox(exi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_aoflux_calc_xao_ox

  !================================================================================================

  function prep_aoflux_get_xao_ox()
    type(mct_aVect), pointer :: prep_aoflux_get_xao_ox(:)
    prep_aoflux_get_xao_ox => xao_ox(:)   
  end function prep_aoflux_get_xao_ox

  function prep_aoflux_get_xao_ax()
    type(mct_aVect), pointer :: prep_aoflux_get_xao_ax(:)
    prep_aoflux_get_xao_ax => xao_ax(:)   
  end function prep_aoflux_get_xao_ax

  function prep_aoflux_get_o2x_ax()
    type(mct_aVect), pointer :: prep_aoflux_get_o2x_ax(:)
    prep_aoflux_get_o2x_ax => o2x_ax(:)   
  end function prep_aoflux_get_o2x_ax

end module prep_aoflux_mod


