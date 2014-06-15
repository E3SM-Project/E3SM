!! This CARMA model is for dust aerosols and is based upon Su & Toon, JGR, 2009;
!! Su & Toon, ACP 2011.
!!
!! These dust are not currently radiatively active and do not replace the dust
!! in CAM; however, this is something that could be done in the future.
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!
!! and adds some local functions used to do sea salt emission:
!!
!!   - CARMA_SurfaceWind()
!!   - CARMA_SurfaceWind()
!!   - WeibullWind()
!!
!! @version July-2012
!! @author  Lin Su, Pengfei Yu, Chuck Bardeen 
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

  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radconstants,   only: nswbands, nlwbands
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
  integer, public, parameter      :: NBIN     = 16              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases

  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_DUST   = 1           !! sea salt composition

  real(kind=f), parameter         :: rClay = 1e-4_f         !! silt/clay particle radius boundary (cm)

  integer                         :: nClay                  !! Number of clay bins (r < 1 um)
  integer                         :: nSilt                  !! Number of silt bins
  real(kind=f)                    :: clay_mf(NBIN)          !! clay mass fraction (fraction)  
  real(kind=f), allocatable, dimension(:,:) :: soil_factor  !! Soil Erosion Factor (fraction) 

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
    integer                            :: LUNOPRT              ! logical unit number for output
    logical                            :: do_print             ! do print output?
    real(kind=f), parameter            :: RHO_DUST = 2.65_f    ! dry density of dust particles (g/cm^3) -Lin Su 
    real(kind=f), parameter            :: rmin     = 1.19e-5_f ! minimum radius (cm)
    real(kind=f), parameter            :: vmrat    = 2.371_f   ! volume ratio
    
    ! Default return code.
    rc = RC_OK    
    
    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_DefineModel: CARMA_Get failed.")
    
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_soilerosion_file = ', carma_soilerosion_file
    end if

    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    call CARMAGROUP_Create(carma, 1, "dust", rmin, vmrat, I_SPHERE, 1._f, .false., &
                          rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="CRDUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, 1, 1, "dust", RHO_DUST, I_INVOLATILE, I_DUST, rc, shortname="CRDUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')
    
    
    ! Define the Solutes
    
    
    ! Define the Gases
    
    
    ! Define the Processes

    
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
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
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
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the bulk mass from the CARMA state.
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version Dec-2010
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use camsrfexch,    only: cam_in_t
    use cam_history,   only: outfld
    
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
    
    integer      :: ilat(pcols)             ! latitude index 
    integer      :: ilon(pcols)             ! longitude index
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: igroup                  ! the index of the carma aerosol group
    character(len=32) :: shortname          ! the shortname of the group
    
    ! -------- local variables added for dust model ------------
    real(r8), parameter :: ch = 0.5e-9_r8                     ! dimensional factor & tuning number as it's model resolution dependent (kgs^2/m^5)!!!
    real(r8)            :: r(NBIN)                            ! bin center (cm)
    real(r8)            :: uth                                ! threshold wind velocity (m/s)

    real(r8)            :: uv10                               ! 10 m wind speed (m/s)
    real(r8)            :: cd10                               ! 10-m drag coefficient ()
    real(r8)            :: wwd                                ! raw wind speed (m/s) 
    real(r8)            :: sp                                 ! mass fraction for soil factor
    integer             :: idustbin                           ! ibin to use for dust production, smallest silt bin for clay
    real(r8)            :: soilfact(pcols)                    ! soil erosion factor (for debug)

    ! Default return code.
    rc = RC_OK

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    call get_lat_all_p(lchnk, ncol, ilat)
    call get_lon_all_p(lchnk, ncol, ilon)

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.            
    tendency(:ncol, :pver) = 0.0_r8
        
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
    if (RC < RC_ERROR) return
    
    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, r=r)
    if (RC < RC_ERROR) return
    
    if (shortname .eq. "CRDUST") then
    
      ! Is this clay or silt?
      !
      ! NOTE: It is assumed that 90% of the mass will be silt and 10% will
      ! be clay.
      !
      ! NOTE: For clay bins, use the smallest silt bin to calculate the
      ! mass and then scale that into each clay bin based upon interpolation of
      ! Tegen and Lacis [1996].
      if (r(ibin) >= rClay) then
        sp         = 0.9_r8 / nSilt
        idustbin   = ibin
      else
        sp         = 0.1_r8 / nClay
        idustbin   = nClay + 1 
      end if

      ! Process each column.
      do icol = 1,ncol
      
        call CARMA_SurfaceWind(carma, state, icol, ilat(icol), ilon(icol), ielem, igroup, idustbin, cam_in, uv10, wwd, uth, rc) 

        ! Is the wind above the threshold for dust production?
        if (uv10 > uth) then
          surfaceFlux(icol) = ch * soil_factor(ilat(icol), ilon(icol)) * sp * &
                              wwd * (uv10 - uth)           
        endif
        
        ! Scale the clay bins based upon the smallest silt bin.   
        surfaceFlux(icol) = clay_mf(ibin) * surfaceFlux(icol)
        
        ! Save off the soil erosion factor so it can be output.
        soilfact(icol) = soil_factor(ilat(icol), ilon(icol))
      end do

      ! For debug purposes, output the soil erosion factor.
      call outfld('CRSLERFC', soilfact, pcols, lchnk)
    end if        
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use cam_history,  only: addfld, add_default, phys_decomp
    use constituents, only: pcnst

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! -------- local variables ----------
    integer            :: ibin                                ! CARMA bin index
    real(r8)           :: r(carma%f_NBIN)                     ! bin center (cm)
    integer            :: count_Silt                          ! count number for Silt
    integer            :: igroup                              ! the index of the carma aerosol group
    integer            :: ielem                               ! the index of the carma aerosol element
    character(len=32)  :: shortname                           ! the shortname of the element
    integer            :: LUNOPRT                             ! logical unit number for output
    logical            :: do_print                            ! do print output?
    
        
    ! Default return code.
    rc = RC_OK

    ! Determine how many clay and how many silt bins there are, based
    ! upon the bin definitions and rClay.
    !
    ! TBD: This should use the radii rather than being hard coded.
    ! nClay = 8
    ! nSilt = NBIN - nClay   
    do ielem = 1, NELEM    
       ! To get particle radius
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
       if (RC < RC_ERROR) return
       
       call CARMAGROUP_GET(carma, igroup, rc, r=r)
       if (RC < RC_ERROR) return
       
       if (shortname .eq. "CRDUST") then
          count_Silt = 0
          do ibin = 1, NBIN
             if (r(ibin) >= rclay) then
                count_Silt = count_Silt + 1
             else
             end if
          end do       
          nSilt = count_Silt
          nClay = NBIN - nSilt     
       end if       
    end do
    
    ! Read in the soil factors.
    call CARMA_ReadSoilErosionFactor(carma, rc)
    if (RC < RC_ERROR) return
    
    ! To determine Clay Mass Fraction
    do ielem = 1, NELEM    
       ! To get particle radius
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
       if (RC < RC_ERROR) return

       if (shortname .eq. "CRDUST") then
          call CARMA_ClayMassFraction(carma, igroup, rc) 
       end if       
    end do
    
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")

      if (do_print) then
        write(carma%f_LUNOPRT,*) 'Initializing CARMA dust model ...'
        write(carma%f_LUNOPRT,*) 'nClay = ', nClay, ' nSilt = ', nSilt
        write(carma%f_LUNOPRT,*) 'clay_mf = ', clay_mf    
        write(carma%f_LUNOPRT,*) 'soil_factor = ', soil_factor
        
        write(carma%f_LUNOPRT,*) 'CARMA dust initialization complete'
      end if
    end if
    
    call addfld('CRSLERFC', 'fraction', 1, 'A', 'CARMA soil erosion factor', phys_decomp)
    
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
    !
    ! NOTE: Initialized to 0. by the caller, so nothing needs to be done.

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
    
    return
  end subroutine CARMA_WetDeposition 


  !!  Determines the mass fraction for the clay (submicron) bins based upon
  !!  Tegen and Lacis [1996]. The total fraction for all clay bins should
  !!  add up to 1.
  !!
  !!  NOTE: WOuld it be better to interpolate this into the bins rather than
  !!  assigning all CARMA bins within a Tegen & Lacis bin the same value?
  !!
  !!  NOTE: Should any mass go to bins smaller than the smallest one used by
  !!  Tegen and Lacis?
  !!
  !!  @version July-2012 
  !!  @author  Lin Su, Pengfei Yu, Chuck Bardeen 
  subroutine CARMA_ClayMassFraction(carma, igroup, rc)
    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: igroup      !! the carma group index
    integer, intent(inout)               :: rc          !! return code, negative indicates failure

    ! Bins and mass fraction from Tegen and Lacis.
    integer, parameter  :: NBIN_TEGEN = 4    
    real(r8)            :: tl_rmin(NBIN_TEGEN) = (/ 1.e-5_r8,  1.8e-5_r8, 3.e-5_r8, 6.e-5_r8 /)
    real(r8)            :: tl_rmax(NBIN_TEGEN) = (/ 1.8e-5_r8, 3.e-5_r8,  6.e-5_r8, 1.e-4_r8 /)
    real(r8)            :: tl_mf(NBIN_TEGEN)   = (/ 0.009_r8,  0.081_r8,  0.234_r8, 0.676_r8 /)

    ! Local Variables
    integer, parameter  :: IBELOW = 1    
    integer, parameter  :: IABOVE = 6    
    integer             :: tl_count(NBIN_TEGEN+2)  ! count number in Tegen and Lacis ranges
    integer             :: ind_up(NBIN_TEGEN+2)
    integer             :: ind_low(NBIN_TEGEN+2)
    integer             :: j                    ! local index number
    integer             :: ibin                 ! carma bin index
    real(r8)            :: r(carma%f_NBIN)      ! CARMA bin center (cm)
     
    ! Default return code.
    rc = RC_OK
    
    ! Interpolate from Tegen and Lacis.
    call CARMAGROUP_GET(carma, igroup, rc, r=r)
    if (RC < RC_ERROR) return
    
    ! Figure out how many of the CARMA bins are in each of the Tegen and Lacis
    ! ranges.
    tl_count(:) = 0
    
    do ibin = 1, NBIN
    
      ! Smaller than the range.
      if (r(ibin) < tl_rmin(1)) then
        tl_count(IBELOW) = tl_count(IBELOW) + 1
      end if
      
      ! In the range
      do j = 1, NBIN_TEGEN
        if (r(ibin) < tl_rmax(j) .and. r(ibin) >= tl_rmin(j)) then
          tl_count(j+1) = tl_count(j+1) + 1
        end if
      end do

      ! Bigger than the range.
      if (r(ibin) >= tl_rmax(NBIN_TEGEN)) then
        tl_count(IABOVE) = tl_count(IABOVE) + 1
      end if       
    end do

    ! Determine where the boundaries are between the TEGEN bins and
    ! the CARMA bin structure.
    ind_up(:)   = 0
    ind_low(:)  = 0
    ind_up (IBELOW)  = tl_count(IBELOW)
    ind_low(IBELOW)  = min(1, tl_count(IBELOW))
    
    do j = 1, 5
      ind_up (j+1) = ind_up(j) + tl_count(j+1)
      ind_low(j+1) = ind_up(j) + min(tl_count(j+1), 1)
    end do
    
    ! No mass to bins smaller than the smallest size.
    clay_mf(:) = 0._r8
    
    ! NOTE: This won't work right if the dust bins are coarser than
    ! the Tegen and Lacis bins. In this case mass fraction would need
    ! to be combined from the Tegen & Lacis bins into a CARMA bin. 
    do j = 1, NBIN_TEGEN
      if (tl_count(j+1) > 0) then
        clay_mf(ind_low(j+1):ind_up(j+1)) = tl_mf(j) / tl_count(j+1)
      end if
    end do 
                                       
    clay_mf(ind_low(IABOVE):) = 1._r8

    return
  end subroutine CARMA_ClayMassFraction

                                              
  !! Calculate the sea surface wind with a Weibull distribution.
  !!
  !! NOTE: This should be combined with a similar routine in the sea salt
  !! model, and any differences should be control by parameters into this
  !! routine (and perhaps namelist variables).
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version July-2012
  subroutine CARMA_SurfaceWind(carma, state, icol, ilat, ilon, ielem, igroup, ibin, cam_in, uv10, wwd, uth, rc)
    use ppgrid,           only: pcols, pver
    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t
 
    implicit none

    ! in and out field
    type(carma_type), intent(in)        :: carma                 !! the carma object
    type(physics_state), intent(in)     :: state                 !! physics state   
    integer, intent(in)                 :: icol                  !! column index
    integer, intent(in)                 :: ilat                  !! latitude index 
    integer, intent(in)                 :: ilon                  !! longitude index
    integer, intent(in)                 :: ielem                 !! element index
    integer, intent(in)                 :: igroup                !! group index
    integer, intent(in)                 :: ibin                  !! bin index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: uv10                  !! the 10m wind speed (m/s)
    real(r8), intent(out)               :: wwd                   !! the 10m wind speed  with Weibull applied (m/s)
    real(r8), intent(out)               :: uth                   !! the 10m wind threshold (m/s)
    integer,  intent(inout)             :: rc                    !! return code, negative indicates failure

    real(r8), parameter                 :: vk = 0.4_r8           ! von Karman constant
    real(r8)                            :: r(NBIN)               ! CARMA bin center (cm)
    real(r8)                            :: rhop(NBIN)            ! CARMA partile element density (g/cm3)
    real(r8)                            :: uthfact               !     
    integer                             :: iepart                ! element in group containing the particle concentration
    real(r8), parameter                 :: rhoa = 1.25e-3_r8     ! Air density at surface
    
    rc = RC_OK
    
    ! Get the 10 meter wind speed
    uv10 = cam_in%u10(icol)

    ! Calculate the threshold wind speed of each bin [Marticorena and Bergametti,1995]
    ! note that in cgs units --> m/s
    call CARMAGROUP_GET(carma, igroup, rc, r=r)
    if (RC < RC_ERROR) return
    
    ! Define particle # concentration element index for current group
    call CARMAELEMENT_Get(carma, ielem, rc, rho=rhop)
    if (RC < RC_ERROR) return
        
    if (cam_in%soilw(icol) > 0._r8 .AND. cam_in%soilw(icol) < 0.5_r8) then
       uthfact = 1.2_r8 + 0.2_r8*log10(cam_in%soilw(icol))
       if (r(ibin) > 2.825e-5_r8) then  ! r(4) = 2.825e-5 cm
           uth = uthfact * 1.e-2_r8 * 0.13_r8 * sqrt(rhop(ibin)*GRAV*r(ibin)*2._r8/rhoa) &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/(r(ibin)*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*(r(ibin)*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       else
           uth = uthfact*1.e-2_r8* 0.13_r8 * sqrt(rhop(ibin)*GRAV*(.75e-4_r8)*2./rhoa)   &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/((.75e-4_r8)*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*((.75e-4_r8)*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       endif
    else
       uth = uv10
    endif    

    ! Use Weibull with Lansing's estimate for shape.
    call WeibullWind(uv10, uth, 2._r8, wwd)

    return
  end subroutine CARMA_SurfaceWind


  !! Read in the dust source (soil) erodibility factor from a NETCDF file. In this
  !! processes, the data is regridded from the source size to the size needed by the
  !! model.
  !!
  !! NOTE: This is currently doing 2-D interpolation, but it really should be doing
  !! regridding.
  !!
  !! @author  Pengfei Yu
  !! @version July-2012
  subroutine CARMA_ReadSoilErosionFactor(carma, rc)
    use pmgrid,        only: plat, plon
    use ioFileMod,     only: getfil
    use wrap_nf
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish    
    
    implicit none

    type(carma_type), intent(in)              :: carma                 !! the carma object
    integer, intent(out)                      :: rc                    !! return code, negative indicates failure

    ! local variables
    integer                                   :: idvar, f_nlon, f_nlat, idlat, idlon
    integer                                   :: fid, fid_lon, fid_lat
    real(r8), allocatable, dimension(:,:)     :: ero_factor, ero_factor1
    character(len=256)                        :: ero_file
    real(r8), allocatable, dimension(:)       :: ero_lat               ! latitude dimension
    real(r8), allocatable, dimension(:)       :: ero_lon               ! latitude dimension
    type (interp_type)                        :: wgt1, wgt2
    real(r8)                                  :: lat(plat), lon(plon)
    integer                                   :: i

    rc = RC_OK

    ! Open the netcdf file (read only)
    call getfil(carma_soilerosion_file, ero_file, 0)
    call wrap_open(ero_file, 0, fid)
  
    ! Get file dimensions
    call wrap_inq_dimid(fid, 'plon', fid_lon)
    call wrap_inq_dimid(fid, 'plat', fid_lat)
    call wrap_inq_dimlen(fid, fid_lon, f_nlon)
    call wrap_inq_dimlen(fid, fid_lat, f_nlat)
  
    allocate(ero_lat(f_nlat))
    allocate(ero_lon(f_nlon))
    allocate(ero_factor (f_nlon, f_nlat))
    allocate(ero_factor1(plon, plat))
    allocate(soil_factor(plat, plon))
    
    ! Read in the tables.
    call wrap_inq_varid(fid, 'new_source', idvar)
    i = nf90_get_var (fid, idvar, ero_factor)
    if (i/=NF90_NOERR) then
       write(iulog,*)'CARMA_ReadSoilErosionFactor: error reading varid =', idvar
       call handle_error (i)
    end if
    call wrap_inq_varid(fid, 'plat', idlat)
    call wrap_get_var_realx(fid, idlat,  ero_lat)
    call wrap_inq_varid(fid, 'plon', idlon)
    call wrap_get_var_realx(fid, idlon,  ero_lon)
            
    ! Close the file.
    call wrap_close(fid)
    
    ! NOTE: Is there a better way to get all of the dimensions
    ! needed for the model grid? Seems like it shouldn't be hard
    ! coded here.
    do i = 1, plat
       lat(i) = 180._r8 / (plat-1) * (i-1) - 90._r8
    end do
    
    do i = 1, plon
       lon(i) = 360._r8 / plon * (i-1)
    end do
    
    call lininterp_init(ero_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(ero_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(ero_factor, f_nlon, f_nlat, ero_factor1, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)
    
    soil_factor(:plat, :plon) = transpose(ero_factor1(:plon, :plat))
    
    deallocate(ero_lat)
    deallocate(ero_lon)
    deallocate(ero_factor)
    deallocate(ero_factor1)
    
    return
  end subroutine CARMA_ReadSoilErosionFactor


  !! Calculate the nth mean of u using Weibull wind distribution
  !! considering the threshold wind velocity. This algorithm
  !! integrates from uth to infinite (u^n P(u)du )
  !!  
  !! @author  Tianyi Fan
  !! @version August-2010
   subroutine WeibullWind(u, uth, n, uwb, wbk)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use shr_spfn_mod, only: gamma =>  shr_spfn_gamma, &
         igamma => shr_spfn_igamma

    implicit none
  
    real(r8), intent(in)  :: u      ! mean wind speed
    real(r8), intent(in)  :: uth    ! threshold velocity
    real(r8), intent(in)  :: n      ! the rank of u in the integration
    real(r8), intent(out) :: uwb    ! the Weibull distribution
    real(r8), intent(in), optional ::  wbk    ! the shape parameter
  
    ! local variable
    real(r8)  :: k                  ! the shape parameter in Weibull distribution
    real(r8)  :: c                  ! the scale parameter in Weibull distribution
  
    if (present(wbk)) then
      k = wbk
    else
      k = 0.94*u**0.5_r8            ! follow Grini and Zender, 2004JGR
 !    k = 2.5_r8                   ! Lansing's estimate
    end if 
  
    ! If u is 0, then k can be 0, which makes a lot of this undefined.
    ! Just return 0. in this case.
    if (u == 0._r8) then
      uwb = 0._r8
    else 
      c   = u * (gamma(1._r8 + 1._r8 / k))**(-1._r8)  
      uwb = c**n * igamma(n / k + 1._r8, (uth / c)**k)
    end if

  end subroutine WeibullWind
  
end module
