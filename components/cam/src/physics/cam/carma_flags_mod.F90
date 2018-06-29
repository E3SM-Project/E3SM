!! This module handles reading the namelist and provides access to some other flags
!! that control CARMA's behavior.
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Aug-2010
module carma_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_readnl                   ! read the carma namelist
  

  ! Namelist flags
  !
  ! NOTE: Setting the carma_flag to false prevents CARMA from doing any microphysics
  ! calculations, but it will still initialize itself. This allows the same build and
  ! namelist to be used, but the CARMA processing diabled. Use the configure option
  ! -carma none to totally disable CARMA and prevent even the register from happening.
  logical, public                :: carma_flag        = .false.   ! If .true. then turn on CARMA microphysics in CAM
  logical, public                :: carma_do_aerosol  = .true.    ! If .true. then CARMA is processed after surface coupling
  logical, public                :: carma_do_cldice   = .false.   ! If .true. then do cloud ice
  logical, public                :: carma_do_cldliq   = .false.   ! If .true. then do cloud liquid
  logical, public                :: carma_do_clearsky = .false.   ! If .true. then do clear sky particle calculations
  logical, public                :: carma_do_coag     = .false.   ! If .true. then do coagulation
  logical, public                :: carma_do_detrain  = .false.   ! If .true. then do detrain
  logical, public                :: carma_do_drydep   = .false.   ! If .true. then do dry deposition
  logical, public                :: carma_do_emission = .false.   ! If .true. then do emission
  logical, public                :: carma_do_fixedinit= .false.   ! If .true. then do fixed initialization to a reference state
  logical, public                :: carma_do_hetchem  = .false.   ! If .true. then CARMA sulfate surface area density used in heterogeneous chemistry
  logical, public                :: carma_do_explised = .false.   ! If .true. then do sedimentation with substepping
  logical, public                :: carma_do_incloud  = .false.   ! If .true. then do incloud particle calculations
  logical, public                :: carma_do_grow     = .false.   ! If .true. then do growth
  logical, public                :: carma_do_optics   = .false.   ! If .true. then do optical properties file
  logical, public                :: carma_do_pheat    = .false.   ! If .true. then do particle heating
  logical, public                :: carma_do_pheatatm = .false.   ! If .true. then do particle heating of atmosphere
  logical, public                :: carma_do_substep  = .false.   ! If .true. then do substeping
  logical, public                :: carma_do_thermo   = .false.   ! If .true. then do solve thermodynamics equation
  logical, public                :: carma_do_wetdep   = .false.   ! If .true. then do wet deposition
  logical, public                :: carma_do_vdiff    = .false.   ! If .true. then do vertical brownian diffusion
  logical, public                :: carma_do_vtran    = .false.   ! If .true. then do vertical transport
  integer, public                :: carma_maxsubsteps = 1         ! Maximum number of time substeps allowed
  integer, public                :: carma_minsubsteps = 1         ! Minimum number of time substeps allowed
  integer, public                :: carma_maxretries  = 8         ! Maximum number of time substeps allowed
  real(r8), public               :: carma_conmax      = 0.1_r8    ! Minumum relative concentration to consider in substep
  real(r8), public               :: carma_dgc_threshold  = 0.0_r8 ! When non-zero, the largest percentage change in gas concentration allowed per substep.
  real(r8), public               :: carma_ds_threshold  = 0.0_r8  ! When non-zero, the largest percentage change in gas saturation allowed per substep.
  real(r8), public               :: carma_dt_threshold  = 0.0_r8  ! When non-zero, the largest change in temperature (K) allowed per substep.
  real(r8), public               :: carma_tstick      = 1.0_r8    ! Thermal accommodation coefficient
  real(r8), public               :: carma_gsticki     = 0.93_r8   ! Growth accommodation coefficient for ice
  real(r8), public               :: carma_gstickl     = 1.0_r8    ! Growth accommodation coefficient for liquid
  real(r8), public               :: carma_cstick      = 1.0_r8    ! Coagulation accommodation coefficient
  real(r8), public               :: carma_rhcrit      = 1.0_r8    ! Critical relative humidity for liquid clouds
  real(r8), public               :: carma_vf_const    = 0.0_r8    ! If specified and non-zero, constant fall velocity for all particles [cm/s]
  character(len=256), public     :: carma_reftfile    = 'carma_reft.nc'  ! path to the file containing the reference temperature profile
  character(len=32), public      :: carma_model       = "none"    ! String (no spaces) that identifies the model

contains


  !! Read the CARMA runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Aug-2010
  subroutine carma_readnl(nlfile)
  
    ! Read carma namelist group.
  
    use cam_abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use carma_model_flags_mod, only: carma_model_readnl
    use perf_mod
  
    ! args
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! local vars
  
    integer :: unitn, ierr
  
    ! read namelist for CARMA
    namelist /carma_nl/ &
      carma_flag, &
      carma_do_aerosol, &
      carma_do_cldliq, &
      carma_do_cldice, &
      carma_do_clearsky, &
      carma_do_coag, &
      carma_do_detrain, &
      carma_do_drydep, &
      carma_do_emission, &
      carma_do_fixedinit, &
      carma_do_hetchem, &
      carma_do_explised, &
      carma_do_incloud, &
      carma_do_grow, &
      carma_do_optics, &
      carma_do_pheat, &
      carma_do_pheatatm, &
      carma_do_substep, &
      carma_do_thermo, &
      carma_do_wetdep, &
      carma_do_vdiff, &
      carma_do_vtran, &
      carma_maxsubsteps, &
      carma_minsubsteps, &
      carma_maxretries, &
      carma_model, &
      carma_reftfile, &
      carma_conmax, &
      carma_dgc_threshold, &
      carma_ds_threshold, &
      carma_dt_threshold, &
      carma_tstick, &
      carma_gsticki, &
      carma_gstickl, &
      carma_cstick, &
      carma_rhcrit, &
      carma_vf_const
  
    logical :: shanlog(22)
    integer :: shanint(3)
    real(r8)  :: shanr8(10)

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
  
#ifdef SPMD
    shanlog(1)  = carma_flag
    shanlog(2)  = carma_do_aerosol
    shanlog(3)  = carma_do_cldliq
    shanlog(4)  = carma_do_cldice
    shanlog(5)  = carma_do_clearsky
    shanlog(6)  = carma_do_coag
    shanlog(7)  = carma_do_detrain
    shanlog(8)  = carma_do_drydep
    shanlog(9)  = carma_do_emission
    shanlog(10) = carma_do_fixedinit
    shanlog(11) = carma_do_hetchem
    shanlog(12) = carma_do_explised
    shanlog(13) = carma_do_incloud
    shanlog(14) = carma_do_grow
    shanlog(15) = carma_do_optics
    shanlog(16) = carma_do_pheat
    shanlog(17) = carma_do_pheatatm
    shanlog(18) = carma_do_substep
    shanlog(19) = carma_do_thermo
    shanlog(20) = carma_do_wetdep
    shanlog(21) = carma_do_vdiff
    shanlog(22) = carma_do_vtran
    call mpibcast (shanlog,            22 ,mpilog, 0,mpicom)
    carma_flag         = shanlog(1)
    carma_do_aerosol   = shanlog(2)
    carma_do_cldliq    = shanlog(3)
    carma_do_cldice    = shanlog(4)
    carma_do_clearsky  = shanlog(5)
    carma_do_coag      = shanlog(6)
    carma_do_detrain   = shanlog(7)
    carma_do_drydep    = shanlog(8)
    carma_do_emission  = shanlog(9)
    carma_do_fixedinit = shanlog(10)
    carma_do_hetchem   = shanlog(11)
    carma_do_explised  = shanlog(12)
    carma_do_incloud   = shanlog(13)
    carma_do_grow      = shanlog(14)
    carma_do_optics    = shanlog(15)
    carma_do_pheat     = shanlog(16)
    carma_do_pheatatm  = shanlog(17)
    carma_do_substep   = shanlog(18)
    carma_do_thermo    = shanlog(19)
    carma_do_wetdep    = shanlog(20)
    carma_do_vdiff     = shanlog(21)
    carma_do_vtran     = shanlog(22)

    shanint(1) = carma_maxsubsteps
    shanint(2) = carma_minsubsteps
    shanint(3) = carma_maxretries
    call mpibcast (shanint,     3 ,mpiint, 0,mpicom)
    carma_maxsubsteps = shanint(1)
    carma_minsubsteps = shanint(2)
    carma_maxretries  = shanint(3)

    shanr8(1) = carma_conmax
    shanr8(2) = carma_dgc_threshold
    shanr8(3) = carma_ds_threshold
    shanr8(4) = carma_dt_threshold
    shanr8(5) = carma_tstick
    shanr8(6) = carma_gsticki
    shanr8(7) = carma_gstickl
    shanr8(8) = carma_cstick
    shanr8(9) = carma_rhcrit
    shanr8(10) = carma_vf_const
    call mpibcast (shanr8,          10 ,mpir8,  0,mpicom)
    carma_conmax         = shanr8(1)
    carma_dgc_threshold  = shanr8(2)
    carma_ds_threshold   = shanr8(3)
    carma_dt_threshold   = shanr8(4)
    carma_tstick         = shanr8(5)
    carma_gsticki        = shanr8(6)
    carma_gstickl        = shanr8(7)
    carma_cstick         = shanr8(8)
    carma_rhcrit         = shanr8(9)
    carma_vf_const       = shanr8(10)
    call mpibcast (carma_model, len(carma_model), mpichar, 0, mpicom)
    call mpibcast (carma_reftfile, len(carma_reftfile), mpichar, 0, mpicom)
#endif

    ! Also cause the CARMA model flags to be read in.
    call carma_model_readnl(nlfile)
  
  end subroutine carma_readnl

end module carma_flags_mod
