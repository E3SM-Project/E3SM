!! This module is used to define a particular CARMA microphysical model. For 
!! simple cases, this may be the only code that needs to be modified. This module
!! defines several constants and has three methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeParticle()
!!
!! These methods define the microphysical model, the particle emissions and
!! the initial conditions of the particles. Each realization of CARMA
!! microphysics has its own version of this file.
!!
!! This file is used to model early earth haze particles. This model is
!! preliminary and used to test the CARMA fractal code. Please talk to
!! Eric Wolf (eric.wolf@colorado.edu) if you are interested in this model.
!!
!! @version May-2013
!! @author  Eric Wolf, Chuck Bardeen 
module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod
  
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use abortutils,     only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition
  
  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 1               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 1               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 40              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases


  !! Relative humidities for mie and radiation calculations. The RRTMG radiation code will interpolate
  !! based upon the current relative humidity from a table built using the specified relative
  !! humidities.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH) = (/ 0._f, 0.5_f, 0.7_f, 0.8_f, 0.9_f, 0.95_f, 0.98_f, 0.99_f /)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_THOLIN       = 1         !! tholin composition

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_THOLIN   = 1         !! tholin aerosol group

  integer, public, parameter      :: I_ELEM_THOLIN  = 1         !! tholin aerosol element
  
  
  ! These variables are all set during initialization and are used to calculate
  ! emission tendencies.
  integer                             :: carma_emis_nLevs        ! number of emission levels
  real(r8), allocatable, dimension(:) :: carma_emis_lev          ! emission levels (Pa)
  real(r8), allocatable, dimension(:) :: carma_emis_rate         ! emission rate lookup table (# cm-3 s-1)
  integer                             :: carma_emis_ilev_min     ! index of minimum level in table 
  integer                             :: carma_emis_ilev_max     ! index of maximum level in table 
  integer                             :: carma_emis_ilev_incr    ! index increment to increase level 
  real(r8)                            :: carma_emis_expected     ! Expected emission rate per column (kg/m2/s)
contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f)                       :: RHO_THOLIN = 0.64     ! density of tholin particles (g/cm)
    real(kind=f), parameter            :: tholin_rmin = 1.e-7_f ! dust minimum radius (cm)
    real(kind=f), parameter            :: tholin_vmrat = 2.5_f  ! dust volume ratio

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?
    real(kind=f)                       :: tholin_rmon = 50e-7_f ! tholin monomer radius (cm)
    
    ! soot fractal dimension
    real(kind=f)                       :: tholin_df(NBIN) = &
           (/ 3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, &
              3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, 3.00000_f, 1.50214_f, 1.50535_f, 1.51331_f, &
              1.53291_f, 1.58003_f, 1.68694_f, 1.89714_f, 2.18998_f, 2.37633_f, 2.39990_f, 2.40000_f, &
              2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, &
              2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f, 2.40000_f /)

    real(kind=f)                       :: tholin_falpha = 1._f  ! soot fractal packing coefficient

    ! Default return code.
    rc = RC_OK

    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.

    call CARMAGROUP_Create(carma, I_GRP_THOLIN, "Tholin", tholin_rmin, tholin_vmrat, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="THOLIN", is_fractal=.true., &
                           rmon=tholin_rmon, df=tholin_df, falpha=tholin_falpha)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')
   
    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_THOLIN, I_GRP_THOLIN, "Tholin", RHO_THOLIN, &
         I_INVOLATILE, I_THOLIN, rc, shortname="THOLIN")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    ! Define the Solutes

    
    ! Define the Gases

    
    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_THOLIN, I_GRP_THOLIN, I_GRP_THOLIN, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    return
  end subroutine CARMA_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
      tnd_qsnow, tnd_nsnow)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none
    
    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    
    ! Default return code.
    rc = RC_OK
        
    return
  end subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)
    use time_manager,     only: is_first_step

    implicit none
    
    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer     :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    real(r8)                             :: mmr(pver) !! elements mass mixing ratio
    integer                              :: ibin      !! bin index
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the mass in each bin from the CAM state.
    
    return
  end subroutine CARMA_DiagnoseBins
  
  
  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)
    
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    real(r8)                             :: mmr(pver) ! mass mixing ration (kg/kg)
    real(r8)                             :: sflx      ! surface flux (kg/m2/s)

    ! Default return code.
    rc = RC_OK

    ! Add the sedimentation and dry deposition fluxes to the hydrophilic black carbon.
    !
    ! NOTE: Don't give the surface model negative values for the surface fluxes.
    ielem = I_ELEM_THOLIN
    do ibin = 1, NBIN
    
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, sedimentationFlux=sflx)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMA_GetBin failed.')
      
      cam_out%ocphidry(icol) = cam_out%ocphidry(icol) + max(sflx, 0._r8)
    end do
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual, is_first_step
    use camsrfexch,    only: cam_in_t
    use tropopause,    only: tropopause_find
    use physconst,     only: gravit
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    integer                            :: ncol                  ! number of columns in chunk
    integer                            :: icol                  ! column index
    integer                            :: igroup                ! the index of the carma aerosol group
    integer                            :: k                     ! vertical index
    integer                            :: ilev                  ! level index in emissions data
    character(len=32)                  :: shortname             ! the shortname of the group
    real(r8)                           :: r(NBIN)               ! bin center
    real(r8)                           :: dr(NBIN)              ! bin width
    real(r8)                           :: rmass(NBIN)           ! bin mass
    real(r8)                           :: pressure              ! pressure (Pa)
    real(r8)                           :: thickness             ! layer thickness (m)
    real(r8)                           :: rate                  ! emission rate (#/cm-3/s)
    real(r8)                           :: massflux              ! emission mass flux (kg/m2/s)
    real(r8)                           :: columnMass            ! mass of the total column (kg/m2/s)
    real(r8)                           :: scale                 ! scaling factor to conserve the expected mass
    
    ! Default return code.
    rc = RC_OK

    ncol = state%ncol

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.            
    tendency(:ncol, :pver) = 0.0_r8


    ! Only do emission for the first bin of the meteor smoke group.
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
    if (RC < RC_ERROR) return
    
    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    ! For meteoritic dust, the source from the smoke only goes into the
    ! smallest bin (~1.3 nm). The depth that the micrometeorite penetrates
    ! is proportional to the pressure, so the emission is a function of
    ! pressure. 
    if ((shortname .eq. "THOLIN") .and. (ibin .eq. 1)) then

      ! Set tendencies for any sources or sinks in the atmosphere.
      do k = 1, pver
        do icol = 1, ncol
      
          pressure = state%pmid(icol, k)
          
          ! This is roughly a log-normal approximation to the production
          ! rate, but only applies from about 70 to 110 km.
          !
          ! NOTE: Based upon US Standard Atmosphere 1976.
          if ((pressure >= carma_emis_lev(carma_emis_ilev_min)) .and. &
              (pressure <= carma_emis_lev(carma_emis_ilev_max))) then

            ! The rates are in terms of # cm-3 s-1, but were really derived
            ! from the mass flux of meteoritic dust. Since we are using a
            ! size different that 1.0 nm for the smallest bin, scale the
            ! number appropriately.
            !
            ! The values are in a lookup table, so find the two numbers
            ! surrounding the pressure and do a linear interpolation on the
            ! rate. This linear search is kind of expensive, particularly if
            ! there are a lot of points.
            ! 
            ! NOTE: The tendency is on a mass mixing ratio (kg/kg/s)
            do ilev = carma_emis_ilev_min, (carma_emis_ilev_max - carma_emis_ilev_incr), carma_emis_ilev_incr
              if ((pressure >= carma_emis_lev(ilev)) .and. (pressure <= carma_emis_lev(ilev+carma_emis_ilev_incr))) then
                rate = carma_emis_rate(ilev)
                
                if (pressure > carma_emis_lev(ilev)) then
                  rate = rate + &
                    ((carma_emis_rate(ilev+carma_emis_ilev_incr) - carma_emis_rate(ilev)) / &
                    (carma_emis_lev(ilev+carma_emis_ilev_incr) - carma_emis_lev(ilev))) * &
                    (pressure - carma_emis_lev(ilev))
                end if
                
                rate = rate * (((1.0e-7_r8)**3) / (r(ibin)**3))
                exit
              end if
            end do
            
            ! Calculate the mass flux in terms of kg/m3/s
            massflux = (rate * rmass(ibin) * 1.0e-3_r8 * 1.0e6_r8)
            
            ! Convert the mass flux to a tendency on the mass mixing ratio.
            thickness = state%zi(icol, k) - state%zi(icol, k+1)
            tendency(icol, k) = (massflux * thickness) / (state%pdel(icol, k) / gravit)        
          end if
        enddo
      enddo

      ! Scale the columns to keep the total mass influx in the column a
      ! constant.
      do icol = 1, ncol
        columnMass = sum(tendency(icol, :) * (state%pdel(icol, :) / gravit))

        ! Protect against divide-by-zero (but not overflow).
        if (columnMass /= 0._r8) then
           scale = carma_emis_expected / columnMass
        else
           scale = 0._r8
        end if

        ! Also apply the relative flux scaling. This needs to be done after
        ! the normalization
        tendency(icol, :) = tendency(icol, :) * scale
      end do
    end if
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use constituents, only: pcnst
    use ioFileMod,    only: getfil
    use wrap_nf
    use mpishorthand

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ilev                  ! level index
    integer                            :: fid                   ! file id
    integer                            :: lev_did               ! level dimension id
    integer                            :: lev_vid               ! level variable id
    integer                            :: rate_vid              ! rate variable
    integer                            :: tmp
    integer                            :: lat_did               ! latitude dimension id
    integer                            :: ltime_did             ! local time dimension id
    integer                            :: time_did              ! time
    integer                            :: lat_vid               ! latitude variable id
    integer                            :: lrf_vid               ! local relative flux variable id
    integer                            :: grf_vid               ! global relative flux variable id
    integer                            :: ltime_vid             ! local time variable id
    character(len=256)                 :: efile                 ! emission file name

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

    ! Default return code.
    rc = RC_OK

    ! Add initialization here.
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
    ! Initialize the emissions rate table.
    if (carma_do_emission) then 
      if (masterproc) then

        ! Open the netcdf file (read only)
        call getfil(carma_emis_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading particle emission rates from ', efile

        call wrap_open(efile, 0, fid)

        ! Alocate the table arrays
        call wrap_inq_dimid(fid, "lev", lev_did)
        call wrap_inq_dimlen(fid, lev_did, carma_emis_nLevs)
      endif
    
#if ( defined SPMD )
      call mpibcast(carma_emis_nLevs, 1, mpiint, 0, mpicom)
#endif

      allocate(carma_emis_lev(carma_emis_nLevs))
      allocate(carma_emis_rate(carma_emis_nLevs))

      if (masterproc) then
        ! Read in the tables.
        call wrap_inq_varid(fid, 'MHAZE', rate_vid)
        call wrap_get_var_realx(fid, rate_vid, carma_emis_rate)

        call wrap_inq_varid(fid, 'lev', lev_vid)
        call wrap_get_var_realx(fid, lev_vid, carma_emis_lev)

        ! Close the file.
        call wrap_close(fid)

        ! Find out where the bounds of the table are and in what order
        ! the pressures levels are in.
        carma_emis_ilev_min = 1
        carma_emis_ilev_max = carma_emis_nLevs

        do ilev = 1, carma_emis_nLevs
          if (carma_emis_rate(ilev) <= 0.0) then
            carma_emis_ilev_min  = ilev + 1
          else
            exit  
          endif
        end do

        do ilev = carma_emis_nLevs, 1, -1
          if (carma_emis_rate(ilev) <= 0.0) then
            carma_emis_ilev_max  = ilev - 1
          else
            exit  
          endif
        end do

        if (carma_emis_lev(carma_emis_ilev_min) < carma_emis_lev(carma_emis_ilev_max)) then
          carma_emis_ilev_incr = 1
        else
          carma_emis_ilev_incr = -1
          tmp = carma_emis_ilev_min
          carma_emis_ilev_min = carma_emis_ilev_max
          carma_emis_iLev_max = tmp 
        endif

        if (do_print) write(LUNOPRT,*) ''
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_nLevs     = ', carma_emis_nLevs
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_min  = ', carma_emis_ilev_min 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_max  = ', carma_emis_ilev_max 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_incr = ', carma_emis_ilev_incr 
        if (do_print) write(LUNOPRT,*) ''
        
        if (do_print) write(LUNOPRT,*) 'level, pressure (Pa), emission rate (# cm-3 sec-1)'
        do ilev = carma_emis_ilev_min, carma_emis_ilev_max, carma_emis_ilev_incr
          if (do_print) write(LUNOPRT,*) ilev, carma_emis_lev(ilev), carma_emis_rate(ilev)
        enddo
        
        if (do_print) write(LUNOPRT, *) 'carma_init(): Total Emission = ', carma_emis_total, ' (kt/yr)'
        carma_emis_expected = ((carma_emis_total * 1e6_r8) / (3600.0_r8 * 24.0_r8 * 365.0_r8)) / &
             (4.0_r8 * PI * ((REARTH / 100._r8) ** 2))
        if (do_print) write(LUNOPRT,*) 'carma_init(): Done with emission table.'

      endif

#if ( defined SPMD )
      call mpibcast(carma_emis_lev,  carma_emis_nLevs, mpir8, 0, mpicom)
      call mpibcast(carma_emis_rate, carma_emis_nLevs, mpir8, 0, mpicom)

      call mpibcast(carma_emis_expected,  1, mpir8,  0, mpicom)

      call mpibcast(carma_emis_ilev_min,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_max,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_incr, 1, mpiint, 0, mpicom)
#endif

    endif
    
    return

    return
  end subroutine CARMA_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeParticle(carma, ielem, ibin, q, gcid, rc)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    real(r8), intent(inout)            :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer, intent(in)                :: gcid(:)  ! global column id
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    
    return
  end subroutine CARMA_InitializeParticle

    
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    ! Add the wet deposition fluxes to the hydrophilic black carbon.
    !
    ! NOTE: Don't give the surface model negative values for the surface fluxes.
    if (ielem == I_ELEM_THOLIN) then
      do icol = 1, state%ncol
        cam_out%ocphiwet(icol) = cam_out%ocphiwet(icol) + max(sflx(icol), 0._r8)
      end do
    end if
    
    return
  end subroutine CARMA_WetDeposition 
  
end module
