module sbetrDriverMod
!
! DESCRIPTION
! module holding subroutines to do out of clm betr application
! created by Jinyun Tang
  use shr_kind_mod , only : r8 => shr_kind_r8
  use BeTR_TimeMod , only : betr_time_type

  implicit none

  private
  save
  public :: sbetrBGC_driver

  character(len=*), parameter :: mod_filename = &
       __FILE__

contains

  subroutine sbetrBGC_driver(base_filename, namelist_buffer)
  !
  !DESCRIPTION
  !driver subroutine for sbetrBGC
  !
  !the rtm is done using the strang splitting approach (Strang, 1968)
  !
  use BeTRSimulationCLM     , only : betr_simulation_clm_type
  use BeTRSimulationStandalone, only : betr_simulation_standalone_type
  use BeTRSimulationALM     , only : betr_simulation_alm_type
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use clm_varpar            , only : nlevtrc_soil
  use decompMod             , only : bounds_type
  use bncdio_pio             , only : file_desc_t
  use clm_instMod           , only : atm2lnd_vars
  use clm_instMod           , only : canopystate_vars
  use clm_instMod           , only : carbonstate_vars, c13state_vars, c14state_vars
  use clm_instMod           , only : carbonflux_vars, c13_cflx_vars, c14_cflx_vars
  use clm_instMod           , only : chemstate_vars, plantMicKinetics_vars
  use clm_instMod           , only : soilhydrology_vars
  use clm_instMod           , only : soilstate_vars
  use clm_instMod           , only : temperature_vars
  use clm_instMod           , only : waterflux_vars
  use clm_instMod           , only : waterstate_vars
  use clm_instMod           , only : cnstate_vars
  use clm_instMod           , only : nitrogenflux_vars, nitrogenstate_vars
  use clm_instMod           , only : phosphorusflux_vars, phosphorusstate_vars
  use clm_instMod           , only : soil_water_retention_curve
  use ColumnType            , only : col
  use clmgridMod            , only : init_clm_vertgrid
  use clm_initializeMod     , only : initialize
  use BeTRSimulation        , only : betr_simulation_type
  use BeTRSimulationFactory , only : create_betr_simulation
  use betr_constants        , only : betr_namelist_buffer_size, betr_string_length_long, betr_filename_length
  use ForcingDataType       , only : ForcingData_type
  use BeTR_GridMod          , only : betr_grid_type
  use TracerParamsMod       , only : tracer_param_init
  use spmdMod               , only : spmd_init
  use LandunitType          , only : lun
  use PatchType             , only : pft
  use landunit_varcon       , only : istsoil
  use clm_varpar            , only : nlevsno, nlevsoi

  implicit none
  !arguments
  character(len=betr_filename_length)      , intent(in) :: base_filename
  character(len=betr_namelist_buffer_size) , intent(in) :: namelist_buffer

  !local variables
  class(betr_simulation_type), pointer :: simulation
  real(r8)                             :: dtime    !model time step
  real(r8)                             :: dtime2   !half of the model time step
  integer                              :: record

  character(len=80)                    :: subname = 'sbetrBGC_driver'

  type(bounds_type)                    :: bounds
  integer                              :: lbj, ubj
  logical :: continue_run
  type(file_desc_t)                      :: ncid
  character(len=betr_string_length_long) :: simulator_name
  character(len=betr_string_length_long) :: restfname
  class(ForcingData_type), allocatable :: forcing_data
  class(betr_grid_type), allocatable :: grid_data
  class(betr_time_type), allocatable :: time_vars

  character(len=256) :: restfile
  integer :: nstep
  logical :: do_standalone=.false.
  logical :: do_alm=.false.
  logical :: do_clm=.false.

  !initialize parameters
  call read_name_list(namelist_buffer, simulator_name, continue_run)
  simulation => create_betr_simulation(simulator_name)

  !set up mask
  bounds%begc = 1
  bounds%endc = 1
  bounds%begp = 1
  bounds%endp = 1
  bounds%begl = 1
  bounds%endl = 1
  bounds%lbj  = 1
  bounds%ubj  = nlevtrc_soil

  !set up grid
  allocate(grid_data)
  call grid_data%Init(namelist_buffer)
  call init_clm_vertgrid(grid_data%nlevgrnd)

  call initialize(bounds)

  allocate(time_vars)
  call time_vars%Init(namelist_buffer)

  !x print*,'read forcing'
  allocate(forcing_data)
  call forcing_data%ReadData(namelist_buffer, grid_data)

  lbj = 1
  ubj = nlevtrc_soil

  lun%itype(1) = istsoil
  col%landunit(1) = 1
  col%gridcell(1) = 1
  col%npfts(1)    = 1
  col%pfti(1)     = 1
  pft%landunit(1) = 1
  pft%column(1)   = 1
  pft%itype(1)    = 2

  call spmd_init
  !print*,'set filters to load initialization data from input'
  call simulation%BeTRSetFilter(maxpft_per_col=1, boffline=.true.)

  if(continue_run)then
    ! print*,'continue from restart file'
    call read_restinfo(restfname, nstep)
    !set current step
    call time_vars%set_nstep(nstep-1)
    ! print*,'back 1 step'
    call time_vars%print_cur_time()
  endif

  !x print*,'obtain waterstate_vars for initilizations that need it'
  call forcing_data%UpdateForcing(grid_data,                                            &
       bounds, lbj, ubj, simulation%num_soilc, simulation%filter_soilc, time_vars, col, &
       pft, atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars    ,   &
       waterflux_vars, temperature_vars, chemstate_vars, simulation%jtops)

  !x print*,'af init update',forcing_data%t_soi(1,:)
  !print*,'initial water state variable output',time_vars%tstep
  call calc_qadv(ubj, simulation%num_soilc, &
       simulation%filter_soilc, waterstate_vars)

  !x print*,'bf sim init'
  !print*,'base_filename:',trim(base_filename)
  call  simulation%Init(bounds, lun, col, pft, waterstate_vars,namelist_buffer, base_filename)
  !x print*,'af sim init'

  select type(simulation)
  class is (betr_simulation_standalone_type)
    print*,'simulation using standalone-betr'
    do_standalone=.true.
  class is (betr_simulation_alm_type)
    print*,'simulation using alm-betr'
    do_alm=.true.
  class is (betr_simulation_clm_type)
    print*,'simulation using clm-betr'
    do_clm=.true.
  class default
  end select

  !read initial condition from restart file is needed
  if(continue_run)then
    call simulation%BeTRRestartOpen(restfname, flag='read', ncid=ncid)
    call simulation%BeTRRestartOffline(bounds, ncid, simulation%num_soilc, simulation%filter_soilc, flag='read')
    call simulation%BeTRRestartClose(ncid)
    !the following aligns forcing data with correct time stamp
    call time_vars%set_nstep(nstep)
    call time_vars%print_cur_time()
    record = 0
  else
    record = -1
    call time_vars%proc_initstep()
  endif

  !x print*,'bf loop'
  do
    record = record + 1
    call simulation%SetClock(dtime=time_vars%get_step_size(), nelapstep=time_vars%get_nstep())
    !x print*,'prepare for diagnosing water flux'
    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, waterstate_vars=waterstate_vars)

    call simulation%PreDiagSoilColWaterFlux(simulation%num_soilc,  simulation%filter_soilc)

    !x print*,'update forcing for betr'
    !set envrionmental forcing by reading foring data: temperature, moisture, atmospheric resistance
    !from either user specified file or clm history file

    call forcing_data%UpdateForcing(grid_data,                                            &
      bounds, lbj, ubj, simulation%num_soilc, simulation%filter_soilc, time_vars, col, pft, &
      atm2lnd_vars, soilhydrology_vars, soilstate_vars,waterstate_vars,                   &
      waterflux_vars, temperature_vars, chemstate_vars, simulation%jtops)

    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, waterstate_vars=waterstate_vars, &
      waterflux_vars=waterflux_vars, soilhydrology_vars = soilhydrology_vars)

    !x print*,'diagnose water flux'
    call simulation%DiagAdvWaterFlux(simulation%num_soilc, &
      simulation%filter_soilc)

    !now assign back waterflux_vars
    call simulation%RetrieveBiogeoFlux(bounds, 1, nlevsoi, waterflux_vars=waterflux_vars)

    !no calculation in the first step
    if(record==0)cycle
    call simulation%BeginMassBalanceCheck(bounds)

    !x print*,'without drainage'
    !the following call could be lsm specific, so that
    !different lsm could use different definitions of input
    !variables, e.g. clm doesn't use cnstate_vars as public variables
    select type(simulation)
    class is (betr_simulation_alm_type)

      call simulation%CalcSmpL(bounds, 1, nlevsoi, simulation%num_soilc, &
        simulation%filter_soilc, temperature_vars%t_soisno_col, &
        soilstate_vars, waterstate_vars, soil_water_retention_curve)

      call simulation%SetBiophysForcing(bounds, col, pft,                               &
        carbonflux_vars=carbonflux_vars,                                                &
        waterstate_vars=waterstate_vars,         waterflux_vars=waterflux_vars,         &
        temperature_vars=temperature_vars,       soilhydrology_vars=soilhydrology_vars, &
        atm2lnd_vars=atm2lnd_vars,               canopystate_vars=canopystate_vars,     &
        chemstate_vars=chemstate_vars,           soilstate_vars=soilstate_vars, &
        cnstate_vars = cnstate_vars)

      call input_substrates((record==1), bounds, col, simulation%num_soilc, simulation%filter_soilc,&
          cnstate_vars, carbonflux_vars,  nitrogenflux_vars, phosphorusflux_vars)

      call simulation%PlantSoilBGCSend(bounds, col, pft, simulation%num_soilc, simulation%filter_soilc,&
        cnstate_vars,  carbonflux_vars, c13_cflx_vars, c14_cflx_vars,  nitrogenflux_vars, phosphorusflux_vars, &
        plantMicKinetics_vars)

    class default
      call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi,               &
        carbonflux_vars=carbonflux_vars,                                                &
        waterstate_vars=waterstate_vars,         waterflux_vars=waterflux_vars,         &
        temperature_vars=temperature_vars,       soilhydrology_vars=soilhydrology_vars, &
        atm2lnd_vars=atm2lnd_vars,               canopystate_vars=canopystate_vars,     &
        chemstate_vars=chemstate_vars,           soilstate_vars=soilstate_vars)
    end select

    call simulation%StepWithoutDrainage(bounds, col, pft)

    select type(simulation)
    class is (betr_simulation_alm_type)
      call simulation%PlantSoilBGCRecv(bounds, col, pft, simulation%num_soilc, simulation%filter_soilc,&
       carbonstate_vars, carbonflux_vars, c13state_vars, c13_cflx_vars, c14state_vars, c14_cflx_vars, &
       nitrogenstate_vars, nitrogenflux_vars, phosphorusstate_vars, phosphorusflux_vars)
    class default
    end select

    !x print*,'with drainge'
    !set forcing variable for drainage
    call simulation%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi,&
       waterflux_vars=waterflux_vars )
    call simulation%StepWithDrainage(bounds, col)

    !x print*,'do mass balance check'
    call simulation%MassBalanceCheck(bounds)

    !specific for water tracer transport
    !call simulation%ConsistencyCheck(bounds, ubj, simulation%num_soilc,    &
    !  simulation%filter_soilc, waterstate_vars)

    !update time stamp
    call time_vars%update_time_stamp()

    !x print*,'write output'
    call simulation%WriteOfflineHistory(bounds, record, simulation%num_soilc,  &
       simulation%filter_soilc, time_vars, waterflux_vars%qflx_adv_col)

    call time_vars%proc_nextstep()

    if(time_vars%its_time_to_write_restart()) then
       !set restfname
       nstep = time_vars%get_nstep()
       write(restfname,'(A,I8.8,A)')trim(base_filename)//'.',nstep,'.rst.nc'
       call write_restinfo(restfname, nstep)

       call simulation%BeTRRestartOpen(restfname, flag='write', ncid=ncid)

       call simulation%BeTRRestartOffline(bounds, ncid, simulation%num_soilc, simulation%filter_soilc, flag='define')

       call simulation%BeTRRestartOffline(bounds, ncid, simulation%num_soilc, simulation%filter_soilc, flag='write')

       call simulation%BeTRRestartClose(ncid)
    endif
    !x print*,'next step'
    if(time_vars%its_time_to_exit()) then
       exit
    end if

  enddo

  call simulation%WriteRegressionOutput(waterflux_vars%qflx_adv_col)
  call forcing_data%destroy()
  deallocate(forcing_data)
end subroutine sbetrBGC_driver

! ----------------------------------------------------------------------

  subroutine read_name_list(namelist_buffer, simulator_name_arg, continue_run)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod        , only : masterproc, mpicom
    use clm_varctl     , only : iulog
    use abortutils     , only : endrun
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size
    use tracer_varcon  , only : advection_on, diffusion_on, reaction_on, ebullition_on, reaction_method
    implicit none
    ! !ARGUMENTS:
    character(len=betr_namelist_buffer_size) , intent(in)  :: namelist_buffer
    character(len=betr_string_length_long)   , intent(out) :: simulator_name_arg
    logical, intent(out) :: continue_run
    !
    ! !LOCAL VARIABLES:
    integer                                :: nml_error
    character(len=*), parameter            :: subname = 'read_name_list'
    character(len=betr_string_length_long) :: simulator_name
    character(len=betr_string_length_long) :: ioerror_msg
    character(len=betr_string_length_long) :: tempstr='daoiga'
    !-----------------------------------------------------------------------

    namelist / sbetr_driver / simulator_name, continue_run

    namelist / betr_parameters /                  &
         reaction_method,                         &
         advection_on, diffusion_on, reaction_on, &
         ebullition_on

    continue_run=.false.
    simulator_name = ''

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=sbetr_driver, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading sbetr_driver namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' sbetr driver :'
       write(stdout, *)
       write(stdout, *) ' sbetr_driver namelist settings :'
       write(stdout, *)
       write(stdout, sbetr_driver)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    reaction_method = ''
    advection_on    = .true.
    diffusion_on    = .true.
    reaction_on     = .true.
    ebullition_on   =.true.

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=betr_parameters, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call endrun(msg="ERROR reading betr_parameters namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr bgc type :'
       write(stdout, *)
       write(stdout, *) ' betr_parameters namelist settings :'
       write(stdout, *)
       write(stdout, betr_parameters)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    simulator_name_arg = simulator_name

  end subroutine read_name_list

  !-------------------------------------------------------------------------------

  subroutine calc_qadv(ubj, numf, filter, waterstate_vars)

    !
    ! description
    ! calculate advective velocity between different layers
    ! this function is now not used, and will be deleted
    ! USES
    use WaterstateType  , only : waterstate_type
    use ColumnType      , only : column_type
    implicit none
    !ARGUMENTS
    integer, intent(in) :: numf
    integer, intent(in) :: filter(:)
    integer, intent(in) :: ubj
    type(waterstate_type)   , intent(inout) :: waterstate_vars

    integer :: j, fc, c

    do j = 1, ubj
       do fc = 1, numf
          c = filter(fc)
          waterstate_vars%h2osoi_liq_old_col(c, j) = waterstate_vars%h2osoi_liq_col(c, j)
          waterstate_vars%h2osoi_ice_old_col(c, j) = waterstate_vars%h2osoi_ice_col(c, j)
       end do
    end do

  end subroutine calc_qadv
  !-------------------------------------------------------------------------------
  subroutine read_restinfo(restfile, nstep)
  implicit none
  character(len=*), intent(out) :: restfile
  integer, intent(out) :: nstep

  integer :: rpt_unit, rpt_error

  rpt_unit=20
  open(unit=rpt_unit,file='rpoint.betr', status='old',&
          action='read', form='formatted', iostat=rpt_error)
  read(rpt_unit,*)restfile,nstep
  print*,'restart file: ', trim(restfile), 'step=',nstep
  close(rpt_unit)

  end subroutine read_restinfo

  !-------------------------------------------------------------------------------
  subroutine write_restinfo(fname, nstep)

  implicit none
  character(len=*), intent(in) :: fname
  integer, intent(in) :: nstep

  integer :: rpt_unit, rpt_error
  character(len=50) :: rpt_filename='rpoint.betr'

  !write rpoint.betr
  rpt_unit=20
  open(unit=rpt_unit,file=trim(rpt_filename), status='replace',&
     action='write', form='formatted', iostat=rpt_error)
  write(rpt_unit,*)trim(fname), nstep
  close(rpt_unit)
  end subroutine write_restinfo

  !-------------------------------------------------------------------------------
  subroutine input_substrates(yesno, bounds, col, num_soilc, filter_soilc,&
    cnstate_vars, carbonflux_vars,  nitrogenflux_vars, phosphorusflux_vars)

  !
  ! DESCRIPTIONS
  ! add substrates for the decompositition test
  use ColumnType     , only : column_type
  use decompMod             , only : bounds_type
  use CNCarbonFluxType, only : carbonflux_type
  use CNNitrogenFluxType, only : nitrogenflux_type
  use PhosphorusFluxType, only : phosphorusflux_type
  use CNStateType, only : cnstate_type
  use clm_varpar, only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use clm_varpar            , only : nlevsoi
  implicit none
  logical           , intent(in)  :: yesno
  type(bounds_type) , intent(in)  :: bounds
  type(column_type) , intent(in)  :: col ! column type
  integer           , intent(in)  :: num_soilc
  integer           , intent(in)  :: filter_soilc(:)
  type(cnstate_type), intent(inout)  :: cnstate_vars
  type(carbonflux_type), intent(inout):: carbonflux_vars
  type(nitrogenflux_type), intent(inout):: nitrogenflux_vars
  type(phosphorusflux_type), intent(inout):: phosphorusflux_vars

  integer :: fc, j, c
  real(r8):: val_c
  real(r8):: val_n, val_p
  real(r8):: tdz(bounds%begc:bounds%endc)

  if(yesno)then
    val_c = 1.e-6_r8
    val_n = 1.e-10_r8
    val_p = 1.e-12_r8
  else
    val_c = 0._r8
    val_n = 0._r8
    val_p = 0._r8
  endif
  do fc = 1, num_soilc
    c = filter_soilc(fc)
    nitrogenflux_vars%ndep_to_sminn_col(c) = val_n
    phosphorusflux_vars%pdep_to_sminp_col(c) = val_p
  enddo
  do j = 1, nlevsoi
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      tdz(c)=sum(col%dz(c,1:nlevsoi))
    enddo
  enddo

  do j = 1, nlevsoi
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      carbonflux_vars%phenology_c_to_litr_met_c_col(c,j) = val_c !gC/m3/s
      carbonflux_vars%phenology_c_to_litr_cel_c_col(c,j) = val_c !gC/m3/s
      carbonflux_vars%phenology_c_to_litr_lig_c_col(c,j) = val_c !gC/m3/s

      nitrogenflux_vars%phenology_n_to_litr_met_n_col(c,j) = val_c/90._r8 !gN/m3/s
      nitrogenflux_vars%phenology_n_to_litr_cel_n_col(c,j) = val_c/90._r8 !gN/m3/s
      nitrogenflux_vars%phenology_n_to_litr_lig_n_col(c,j) = val_c/90._r8 !gN/m3/s

      phosphorusflux_vars%phenology_p_to_litr_met_p_col(c,j) = val_c/1600._r8 !gP/m3/s
      phosphorusflux_vars%phenology_p_to_litr_cel_p_col(c,j) = val_c/1600._r8 !gP/m3/s
      phosphorusflux_vars%phenology_p_to_litr_lig_p_col(c,j) = val_c/1600._r8 !gP/m3/s

      cnstate_vars%ndep_prof_col(c,j) = col%dz(c,j)/tdz(c)
      cnstate_vars%pdep_prof_col(c,j) = col%dz(c,j)/tdz(c)

    enddo
  enddo

  end subroutine input_substrates
  !-------------------------------------------------------------------------------



end module sbetrDriverMod
