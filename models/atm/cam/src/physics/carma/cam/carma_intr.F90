!! This module is a coupler between the CAM model and the Community Aerosol
!! and Radiation Model for Atmospheres (CARMA) microphysics model. It adds the
!! capabilities of CARMA to CAM, allowing for binned microphysics studies
!! within the CAM framework. This module supports the CAM physics interface, and
!! uses the CARMA and CARMASTATE objects to perform the microphysical
!! calculations.
!!
!! @author  Chuck Bardeen
!! @version July 2009
module carma_intr

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carma_flags_mod
  use carma_model_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use pmgrid,         only: plat, plev, plevp, plon
  use ppgrid,         only: pcols, pver, pverp
  use ref_pres,       only: pref_mid, pref_edge, pref_mid_norm, psurf_ref
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init, &
                            set_dry_to_wet, physics_state_copy
  use phys_grid,      only: get_lat_all_p
  use physconst,      only: avogad, cpair
  use constituents,   only: pcnst, cnst_add, cnst_get_ind, &
                            cnst_name, cnst_longname, cnst_type
  use chem_surfvals,  only: chem_surfvals_get
  use abortutils,     only: endrun
  use physics_buffer, only: physics_buffer_desc, pbuf_add_field, pbuf_old_tim_idx, &
                            pbuf_get_index, pbuf_get_field, dtype_r8


#if ( defined SPMD )
  use mpishorthand
#endif  
  
  implicit none
  
  private
  save

  ! Public interfaces
  
  ! CAM Physics Interface
  public carma_register                 ! register consituents
  public carma_is_active                ! retrns true if this package is active (microphysics = .true.)
  public carma_implements_cnst          ! returns true if consituent is implemented by this package
  public carma_init_cnst                ! initialize constituent mixing ratios, if not read from initial file
  public carma_init                     ! initialize timestep independent variables
  public carma_final                    ! finalize the CARMA module
  public carma_timestep_init            ! initialize timestep dependent variables
  public carma_timestep_tend            ! interface to tendency computation
  public carma_accumulate_stats         ! collect stats from all MPI tasks
  
  ! Other Microphysics
  public carma_emission_tend            ! calculate tendency from emission source function
  public carma_wetdep_tend              ! calculate tendency from wet deposition
  

  ! Private data
  
  ! Particle Group Statistics
  
  ! Gridbox average
  integer, parameter             :: NGPDIAGS         = 12         ! Number of particle diagnostics ...
  integer, parameter             :: GPDIAGS_ND       = 1          ! Number density
  integer, parameter             :: GPDIAGS_AD       = 2          ! Surface area density
  integer, parameter             :: GPDIAGS_MD       = 3          ! Mass density
  integer, parameter             :: GPDIAGS_RE       = 4          ! Effective Radius
  integer, parameter             :: GPDIAGS_RM       = 5          ! Mitchell [2002] Effective Radius
  integer, parameter             :: GPDIAGS_JN       = 6          ! Nucleation Rate
  integer, parameter             :: GPDIAGS_MR       = 7          ! Mass Mixing Ratio
  integer, parameter             :: GPDIAGS_EX       = 8          ! Extinction
  integer, parameter             :: GPDIAGS_OD       = 9          ! Optical Depth
  integer, parameter             :: GPDIAGS_VM       = 10         ! Mass Weighted Fall Velocity
  integer, parameter             :: GPDIAGS_PA       = 11         ! Projected Area
  integer, parameter             :: GPDIAGS_AR       = 12         ! Area Ratio

  ! Particle Bin (Element) Statistics
  integer, parameter             :: NBNDIAGS         = 1          ! Number of bin surface diagnostics ...
  integer, parameter             :: BNDIAGS_TP       = 1          ! Delta Particle Temperature [K]
  
  ! Surface
  integer, parameter             :: NSBDIAGS         = 2          ! Number of bin surface diagnostics ...
  integer, parameter             :: SBDIAGS_DD       = 1          ! Dry deposition flux [kg/m2/s]
  integer, parameter             :: SBDIAGS_VD       = 2          ! Dry deposition velocity [cm/s]
  
  
  ! Gas Statistics
  integer, parameter             :: NGSDIAGS         = 5          ! Number of gas diagnostics ...
  integer, parameter             :: GSDIAGS_SI       = 1          ! saturation wrt ice
  integer, parameter             :: GSDIAGS_SL       = 2          ! saturation wrt water
  integer, parameter             :: GSDIAGS_EI       = 3          ! equilibrium vp wrt ice
  integer, parameter             :: GSDIAGS_EL       = 4          ! equilibrium vp wrt water
  integer, parameter             :: GSDIAGS_WT       = 5          ! weight percent composition for aerosols
  
  ! Step Statistics
  integer, parameter             :: NSPDIAGS         = 2          ! Number of step diagnostics ...
  integer, parameter             :: SPDIAGS_NSTEP    = 1          ! number of substeps
  integer, parameter             :: SPDIAGS_LNSTEP   = 2          ! ln(number of substeps)
  
  ! Defaults not in the namelist
  character(len=10), parameter   :: carma_mixtype     = 'wet'     ! mixing ratio type for CARMA constituents
  integer                        :: LUNOPRT           = -1        ! lun for output
  
  ! Constituent Mappings      
  integer                        :: icnst4elem(NELEM, NBIN)       ! constituent index for a carma element
  integer                        :: icnst4gas(NGAS)               ! constituent index for a carma gas

  character(len=16)              :: btndname(NGROUP, NBIN)        ! names of group per bin tendencies
  character(len=16)              :: etndname(NELEM, NBIN)         ! names of element tendencies
  character(len=16)              :: gtndname(NGAS)                ! names of gas tendencies
  
  ! Flags to indicate whether each constituent could have a CARMA tendency.
  logical                        :: lq_carma(pcnst)
  
  ! The CARMA object stores the configuration inforamtion about CARMA, only one is
  ! is needed per MPI task. In the future, this could potentially be turned into one
  ! per model to allow multiple models with different numbers of bins, ... to be
  ! run simultaneously. However, it is more complicated than that, since some of the
  ! globals here would need to be put into the carma or another object and some sort
  ! of callback mechanism is needed to call the correct model implementations for
  ! the model specific methods.
  type(carma_type), target       :: carma                         ! the carma object


  ! Physics Buffer Indicies  
  integer                        :: ipbuf4gas(NGAS)               ! physics buffer index for a carma gas
  integer                        :: ipbuf4t                       ! physics buffer index for a carma temperature
  integer                        :: ipbuf4sati(NGAS)              ! physics buffer index for a carma saturation over ice
  integer                        :: ipbuf4satl(NGAS)              ! physics buffer index for a carma saturation over liquid
  
  ! Globals used for a reference atmosphere.
  real(kind=f)                   :: carma_t_ref(pver)             ! midpoint temperature (Pa)
  real(kind=f)                   :: carma_h2o_ref(pver)           ! h2o mmmr (kg/kg)
  real(kind=f)                   :: carma_h2so4_ref(pver)         ! h2so4 mmr (kg/kg)


  ! Globals used for total statistics
  real(kind=f)          :: glob_max_nsubstep = 0._f
  real(kind=f)          :: glob_max_nretry   = 0._f
  real(kind=f)          :: glob_nstep        = 0._f
  real(kind=f)          :: glob_nsubstep     = 0._f
  real(kind=f)          :: glob_nretry       = 0._f

  real(kind=f)          :: step_max_nsubstep = 0._f
  real(kind=f)          :: step_max_nretry   = 0._f
  real(kind=f)          :: step_nstep        = 0._f
  real(kind=f)          :: step_nsubstep     = 0._f
  real(kind=f)          :: step_nretry       = 0._f


contains


  !! Read the names of the constituents from CARMA and automatically create
  !! a list of names based on the constituents and the number of size
  !! bins. A naming convention is used to map from CARMA constiuent & bin to
  !! CAM constituent, with the smallest bin being <shortname>01, then next
  !! shortname<02>, ...
  !!
  !! A check is done to see if the CARMA gases are already present. If so,
  !! they gases are linked; otherwise, the gas is added to CAM. The shortname
  !! of the gas is used as the constituent name.
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author Chuck Bardeen
  !! @version May-2009
  subroutine carma_register
    use radconstants,    only : nswbands, nlwbands, &
         get_sw_spectral_boundaries, get_lw_spectral_boundaries
    use cam_logfile,     only : iulog
    use cam_control_mod, only : nsrest  
    use physconst,    only: gravit, p_rearth=>rearth, mwdry, mwh2o
    use phys_control, only: phys_getopts

    implicit none

    integer           :: ielem                    ! element index
    integer           :: ibin                     ! bin index
    integer           :: igas                     ! gas index
    integer           :: igroup                   ! group index
    integer           :: rc                       ! CARMA return code
    character(len=8)  :: c_name                   ! constituent name
    character(len=50) :: c_longname               ! constituent long name
    real(r8)          :: wave(NWAVE)              ! wavelength band centers (cm)
    real(r8)          :: dwave(NWAVE)             ! wavelength band width (cm)
    logical           :: do_wave_emit(NWAVE)      ! do emission in band?
    real(r8)          :: r(NBIN)                  ! particle radius (cm)
    real(r8)          :: rmass(NBIN)              ! particle mass (g)
    character(len=8)  :: shortname                ! short (CAM) name
    character(len=8)  :: grp_short                ! group short (CAM) name
    character(len=50) :: name                     ! long (CARMA) name
    real(r8)          :: wtmol                    ! gas molecular weight
    integer           :: cnsttype                 ! constituent type
    integer           :: maxbin                   ! last prognostic bin

    character(len=16) :: radiation_scheme         ! CAM's radiation package.

    ! Initialize the return code.
    rc = 0
    
    ! Some constants are set on the fly in CAM, so initialize them and any derived "constants" here.
    ! Some of them are needed in CARMA_DefineModel and CARMA_Initialize.
    GRAV = gravit * RM2CGS
    REARTH  = p_rearth * RM2CGS 
    WTMOL_AIR = mwdry 
    WTMOL_H2O = mwh2o 
    R_AIR = RGAS / WTMOL_AIR
    RKAPPA = R_AIR / CP

    ! Setup the lun for output.
    LUNOPRT = iulog

    ! Find out which radiation scheme is active.
    call phys_getopts(radiation_scheme_out = radiation_scheme)
   
    ! Get the wavelength centers for the CAM longwave and shortwave bands
    ! from the radiation code.

    ! Can do this 'in place'; set wave to lower boundaries of all bands,
    ! dwave to upper boundaries.
    call get_lw_spectral_boundaries( wave(:nlwbands   ), dwave(:nlwbands   ), 'cm')
    call get_sw_spectral_boundaries( wave( nlwbands+1:), dwave( nlwbands+1:), 'cm')

    ! Now make dwave the difference and wave the average
    dwave = dwave - wave
    wave  = wave + (dwave / 2._r8)

    ! NOTE: RRTMG does not calculate emission in the shortwave bands and the first and
    ! last shortwave bands overlap the longwave bands. At least the first and last bands
    ! needs to be excluded and perhaps all of the shortwave bands need to be excluded or
    ! double counting will happen for particle emission. The details of this should
    ! probably be moved into the radiation code.
    do_wave_emit(:nlwbands) = .TRUE.
!    do_wave_emit(nlwbands+1:) = .FALSE.
    do_wave_emit(nlwbands+1:) = .TRUE.
    do_wave_emit(nlwbands+1)  = .FALSE.
    do_wave_emit(NWAVE)       = .FALSE.

    ! Create the CARMA object that will contain all the information about the
    ! how CARMA is configured.
      
    call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
         LUNOPRT=LUNOPRT, wave=wave, dwave=dwave, do_wave_emit=do_wave_emit)
    if (rc < 0) call endrun('carma_register::CARMA_Create failed.')
    
    ! Define the microphysical model.
    call CARMA_DefineModel(carma, rc)
    if (rc < 0) call endrun('carma_register::CARMA_DefineModel failed.')
    
    if (masterproc) then
      write(LUNOPRT,*) ''
      write(LUNOPRT,*) 'CARMA general settings for ', trim(carma_model), ' model : '
      write(LUNOPRT,*) '  carma_do_aerosol    = ', carma_do_aerosol
      write(LUNOPRT,*) '  carma_do_cldice     = ', carma_do_cldice
      write(LUNOPRT,*) '  carma_do_cldliq     = ', carma_do_cldliq
      write(LUNOPRT,*) '  carma_do_clearsky   = ', carma_do_clearsky
      write(LUNOPRT,*) '  carma_do_coag       = ', carma_do_coag
      write(LUNOPRT,*) '  carma_do_detrain    = ', carma_do_detrain
      write(LUNOPRT,*) '  carma_do_drydep     = ', carma_do_drydep
      write(LUNOPRT,*) '  carma_do_emission   = ', carma_do_emission
      write(LUNOPRT,*) '  carma_do_fixedinit  = ', carma_do_fixedinit
      write(LUNOPRT,*) '  carma_do_grow       = ', carma_do_grow
      write(LUNOPRT,*) '  carma_do_explised   = ', carma_do_explised
      write(LUNOPRT,*) '  carma_do_incloud    = ', carma_do_incloud
      write(LUNOPRT,*) '  carma_do_pheat      = ', carma_do_pheat
      write(LUNOPRT,*) '  carma_do_pheatatm   = ', carma_do_pheatatm
      write(LUNOPRT,*) '  carma_do_substep    = ', carma_do_substep
      write(LUNOPRT,*) '  carma_do_thermo     = ', carma_do_thermo
      write(LUNOPRT,*) '  carma_do_vdiff      = ', carma_do_vdiff
      write(LUNOPRT,*) '  carma_do_vtran      = ', carma_do_vtran
      write(LUNOPRT,*) '  carma_do_wetdep     = ', carma_do_wetdep
      write(LUNOPRT,*) '  carma_dgc_threshold = ', carma_dgc_threshold
      write(LUNOPRT,*) '  carma_ds_threshold  = ', carma_ds_threshold
      write(LUNOPRT,*) '  carma_dt_threshold  = ', carma_dt_threshold
      write(LUNOPRT,*) '  carma_cstick        = ', carma_cstick
      write(LUNOPRT,*) '  carma_gsticki       = ', carma_gsticki
      write(LUNOPRT,*) '  carma_gstickl       = ', carma_gstickl
      write(LUNOPRT,*) '  carma_tstick        = ', carma_tstick
      write(LUNOPRT,*) '  carma_rhcrit        = ', carma_rhcrit
      write(LUNOPRT,*) '  carma_conmax        = ', carma_conmax
      write(LUNOPRT,*) '  carma_minsubsteps   = ', carma_minsubsteps
      write(LUNOPRT,*) '  carma_maxsubsteps   = ', carma_maxsubsteps
      write(LUNOPRT,*) '  carma_maxretries    = ', carma_maxretries
      write(LUNOPRT,*) '  carma_vf_const      = ', carma_vf_const
      write(LUNOPRT,*) '  cldfrc_incloud      = ', CLDFRC_INCLOUD
      write(LUNOPRT,*) '  carma_reftfile      = ', trim(carma_reftfile)
      write(LUNOPRT,*) ''
    endif
    
    ! Intialize the model based upon the namelist configuration.
    !
    ! NOTE: When used with CAM, the latents heats (of melting and evaporation)
    ! need to be constant for energy balance to work. This allows them to match the
    ! assumptions made in the CAM energy checking and microphysics code.
    call CARMA_Initialize(carma, &
                          rc, &
                          do_clearsky   = carma_do_clearsky, &
                          do_cnst_rlh   = .true., &
                          do_coag       = carma_do_coag, &
                          do_detrain    = carma_do_detrain, &
                          do_drydep     = carma_do_drydep, &
                          do_fixedinit  = carma_do_fixedinit, &
                          do_grow       = carma_do_grow, &
                          do_explised   = carma_do_explised, &
                          do_incloud    = carma_do_incloud, &
                          do_pheat      = carma_do_pheat, &
                          do_pheatatm   = carma_do_pheatatm, &
                          do_print_init = masterproc, &
                          do_substep    = carma_do_substep, &
                          do_thermo     = carma_do_thermo, &
                          do_vdiff      = carma_do_vdiff, &
                          do_vtran      = carma_do_vtran, &
                          minsubsteps   = carma_minsubsteps, &
                          maxsubsteps   = carma_maxsubsteps, &
                          maxretries    = carma_maxretries, &
                          vf_const      = carma_vf_const, &
                          conmax        = carma_conmax, &
                          cstick        = carma_cstick, &
                          dt_threshold  = carma_dt_threshold, &
                          gsticki       = carma_gsticki, &
                          gstickl       = carma_gstickl, &
                          tstick        = carma_tstick)
    if (rc < 0) call endrun('carma_register::CARMA_Initialize failed.')
    
    
    ! The elements and gases from CARMA need to be added as constituents in
    ! CAM (if they don't already exist). For the elements, each radius bin
    ! needs to be its own constiuent in CAM.
    !
    ! Some rules about the constituents:
    !   1) The shortname must be 8 characters or less, so the CARMA short name
    !      is limited to 6 characters and 2 characters for the bin number.
    !   2) The molecular weight is in kg/kmol.
    !   3) The specific heat at constant pressure is in J/kg/K.
    !   4) The consituents are added sequentially.
    
    ! Add a CAM constituents for each bin of each element.
    do ielem = 1, NELEM
    
      call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, shortname=shortname, name=name)
      if (rc < 0) call endrun('carma_register::CARMAELEMENT_Get failed.')
      
      call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, r=r, rmass=rmass, maxbin=maxbin, shortname=grp_short)
      if (rc < 0) call endrun('carma_register::CARMAGROUP_Get failed.')
      
      ! For prognostic groups, all of the bins need to be represented as actual CAM
      ! constituents. Diagnostic groups are determined from state information that
      ! is already present in CAM, and thus their bins only exist in CARMA.
      if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then
        
        do ibin = 1, NBIN
        
          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the paerent model.
          if (ibin <= maxbin) then
    
            write(btndname(igroup, ibin), '(A, I2.2)') trim(grp_short), ibin

            write(c_name, '(A, I2.2)') trim(shortname), ibin
            write(c_longname, '(A, e11.4, A)') trim(name) // ', ', r(ibin)*1.e4_r8, ' um'
             
            ! The molecular weight seems to be used for molecular diffusion, which
            ! doesn't make sense for particles. The CAM solvers are unstable if the 
            ! mass provided is large.
            call cnst_add(c_name, WTMOL_AIR, cpair, 0._r8, icnst4elem(ielem, ibin), &
              longname=c_longname, mixtype=carma_mixtype, is_convtran1=is_convtran1(igroup))
          end if
        end do
      end if
    end do
     
    ! Find the constituent for the gas or add it if not found.
    do igas = 1, NGAS
  
      call CARMAGAS_Get(carma, igas, rc, shortname=shortname, name=name, wtmol=wtmol)
      if (rc < 0) call endrun('carma_register::CARMAGAS_Get failed.')
        
      ! Is the gas already defined?
      call cnst_get_ind(shortname, icnst4gas(igas))
  
      ! For substepping, we need to store the last mmr values for the gas.
      call pbuf_add_field('CG' // shortname, 'global',dtype_r8, (/pcols, pver/), ipbuf4gas(igas))
  
      ! For substepping, we need to store the last supersaturations.
      call pbuf_add_field('CI' // shortname, 'global',dtype_r8, (/pcols, pver/), ipbuf4sati(igas))
      call pbuf_add_field('CL' // shortname, 'global',dtype_r8, (/pcols, pver/), ipbuf4satl(igas))
    end do
     
  
    ! For substepping, we need to store the temperature.
    call pbuf_add_field('CT', 'global',dtype_r8, (/pcols, pver/), ipbuf4t)
    

    ! Create the optical properties files needed for RRTMG radiative transfer
    ! calculations.
    !
    ! NOTE: This only needs to be done once at the start of the run and does not need
    ! to be done for restarts.
    !
    ! NOTE: We only want to do this with RRTMG. If CAM_RT is being used, then skip this.
    if ((masterproc) .and. (nsrest == 0) .and. (radiation_scheme == "rrtmg") .and. (carma_do_optics)) then
      call CARMA_CreateOpticsFile(carma, rc)
       if (rc < 0) call endrun('carma_register::carma_CreateOpticsFiles failed.')
    end if
    
    return
  end subroutine carma_register


  !! Returns true if the CARMA package is active
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author  Chuck Bardeen
  !! @version May 2009
  function carma_is_active()
    implicit none
  
    logical :: carma_is_active
  
    carma_is_active = carma_flag
    
    return
  end function carma_is_active


  !!  Returns true if specified constituent is implemented by CARMA
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author  Chuck Bardeen
  !! @version May 2009
  function carma_implements_cnst(name)
    implicit none
    
    character(len=*), intent(in) :: name   !! constituent name
    logical :: carma_implements_cnst       ! return value
    
    integer :: igroup
    integer :: ielem
    integer :: ibin
    integer :: igas
    integer :: rc
    
    integer :: cnsttype     ! constituent type
    integer :: maxbin       ! last prognostic bin

    rc = 0
    
    carma_implements_cnst = .false.
    
    ! Check each bin to see if it this constituent.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMAELEMENT_Get(carma,  ielem, rc, igroup=igroup)
        if (rc < 0) call endrun('carma_init::CARMAELEMENT_Get failed.')
        
        call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, maxbin=maxbin)
        if (rc < 0) call endrun('carma_init::CARMAGROUP_Get failed.')
        
        if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then
        
          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the parent model.
          if (ibin <= maxbin) then
    
            if (name == cnst_name(icnst4elem(ielem, ibin))) then
              carma_implements_cnst = .true.
              return
            end if
          end if
        end if
      end do
    end do 
    
    ! Check each gas to see if it this constituent.
    do igas = 1, NGAS
      if (name == cnst_name(icnst4gas(igas))) then
         carma_implements_cnst = .true.
         return
      end if
    end do
    
    return
  end function carma_implements_cnst
  

  !! Initialize items in CARMA that only need to be initialized once. This
  !! routine is called after carma_register has been called.
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author  Chuck Bardeen
  !! @version May 2009
  subroutine carma_init
    use cam_history,  only: addfld, add_default, phys_decomp
    use ioFileMod,    only : getfil
    use wrap_nf
    use time_manager, only: is_first_step
  
    implicit none
    
    integer           :: iz           ! vertical index
    integer           :: ielem        ! element index
    integer           :: ibin         ! bin index
    integer           :: igas         ! gas index
    integer           :: igroup       ! group index
    integer           :: icnst        ! constituent index
    integer           :: rc           ! CARMA return code
    character(len=8)  :: sname        ! short (CAM) name
    integer           :: cnsttype     ! constituent type
    integer           :: maxbin       ! last prognostic bin
    logical           :: is_cloud     ! is the group a cloud?
    logical           :: do_drydep    ! is dry deposition enabled?
 
    integer                    :: i
    integer                    :: ier
    integer                    :: ncid, dimid_lev, lev, vid_T
    logical                    :: lexist
    character(len=256)         :: locfn
    integer                    :: nlev
    integer                    :: LUNOPRT              ! logical unit number for output
    logical                    :: do_print             ! do print output?
    

1   format(a6,4x,a11,4x,a11,4x,a11)
2   format(i6,4x,3(1pe11.3,4x))

    ! Initialize the return code.
    rc = 0
    
    ! Set names of constituent sources and declare them as history variables; howver,
    ! only prognostic variables have.
    lq_carma(:) = .false.
     
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMAELEMENT_Get(carma,  ielem, rc, igroup=igroup)
        if (rc < 0) call endrun('carma_init::CARMAELEMENT_Get failed.')
        
        call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, maxbin=maxbin, do_drydep=do_drydep)
        if (rc < 0) call endrun('carma_init::CARMAGROUP_Get failed.')
        
        if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then

          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the parent model.
          if (ibin <= maxbin) then
    
            icnst = icnst4elem(ielem, ibin)
            
            ! Indicate that this is a constituent whose tendency could be changed by
            ! CARMA.
            lq_carma(icnst) = .true.
            
            etndname(ielem, ibin) = trim(cnst_name(icnst))
            
            call addfld(cnst_name(icnst),             'kg/kg   ', pver, 'A', cnst_longname(icnst), phys_decomp)
            call add_default(cnst_name(icnst), 1, ' ')
    
            call addfld(trim(etndname(ielem, ibin))//'TC', 'kg/kg/s ', pver, 'A', &
                 trim(cnst_name(icnst)) // ' tendency', phys_decomp)
            call addfld(trim(etndname(ielem, ibin))//'SF', 'kg/m2/s ', 1,    'A', &
                 trim(cnst_name(icnst)) // ' surface emission', phys_decomp)
            call addfld(trim(etndname(ielem, ibin))//'EM', 'kg/kg/s ', pver, 'A', &
                 trim(cnst_name(icnst)) // ' emission', phys_decomp)
            call addfld(trim(etndname(ielem, ibin))//'WD', 'kg/kg/s ', pver, 'A', &
                 trim(cnst_name(icnst)) // ' wet deposition', phys_decomp)
            call addfld(trim(etndname(ielem, ibin))//'SW', 'kg/m2/s ', 1,    'A', &
                 trim(cnst_name(icnst)) // ' wet deposition flux at surface', phys_decomp)

            if (do_drydep) then
              call addfld(trim(etndname(ielem, ibin))//'DD', 'kg/m2/s ', 1,  'A', &
                   trim(cnst_name(icnst)) // ' dry deposition', phys_decomp)
            end if

            if (carma_do_pheat) then
              call addfld(trim(etndname(ielem, ibin))//'TP', 'K     ', pver, 'A', &
                   trim(cnst_name(icnst)) // ' delta particle temperature', phys_decomp)
            end if
          end if
        end if
      end do
    end do

    do igroup = 1, NGROUP
      call CARMAGROUP_Get(carma, igroup, rc, shortname=sname, is_cloud=is_cloud, do_drydep=do_drydep)
      if (rc < 0) call endrun('carma_init::CARMAGROUP_GetGroup failed.')
      
      ! Gridbox average
      !
      ! NOTE: Would like use flag_xf_fill for the reffective radius fields, but cam_history
      ! currently only supports fill values in the entire column.
      call addfld(trim(sname)//'ND', '#/cm3   ', pver, 'A', trim(sname) // ' number density', phys_decomp)
      call addfld(trim(sname)//'AD', 'um2/cm3 ', pver, 'A', trim(sname) // ' surface area density', phys_decomp)
      call addfld(trim(sname)//'MD', 'g/cm3   ', pver, 'A', trim(sname) // ' mass density', phys_decomp)
      call addfld(trim(sname)//'RE', 'um      ', pver, 'A', trim(sname) // ' effective radius', phys_decomp)
      call addfld(trim(sname)//'RM', 'um      ', pver, 'A', trim(sname) // ' Mitchell effective radius', phys_decomp)
      call addfld(trim(sname)//'JN', '#/cm3/s ', pver, 'A', trim(sname) // ' nucleation rate', phys_decomp)
      call addfld(trim(sname)//'MR', 'kg/kg   ', pver, 'A', trim(sname) // ' mass mixing ratio', phys_decomp)
      call addfld(trim(sname)//'EX', 'km-1    ', pver, 'A', trim(sname) // ' extinction', phys_decomp)
      call addfld(trim(sname)//'OD', '        ', pver, 'A', trim(sname) // ' optical depth', phys_decomp)
      call addfld(trim(sname)//'PA', 'cm2     ', pver, 'A', trim(sname) // ' projected area', phys_decomp)
      call addfld(trim(sname)//'AR', '        ', pver, 'A', trim(sname) // ' area ratio', phys_decomp)
      call addfld(trim(sname)//'VM', 'm/s     ', pver, 'A', trim(sname) // ' fall velocity', phys_decomp)

      call add_default(trim(sname)//'ND', 1, ' ')
      call add_default(trim(sname)//'AD', 1, ' ')
      call add_default(trim(sname)//'MD', 1, ' ')
      call add_default(trim(sname)//'RE', 1, ' ')
      call add_default(trim(sname)//'RM', 1, ' ')
      call add_default(trim(sname)//'MR', 1, ' ')
      call add_default(trim(sname)//'EX', 1, ' ')
      call add_default(trim(sname)//'OD', 1, ' ')
      call add_default(trim(sname)//'PA', 1, ' ')
      call add_default(trim(sname)//'AR', 1, ' ')
      call add_default(trim(sname)//'VM', 1, ' ')
      
      if (carma_do_grow) then
        call add_default(trim(sname)//'JN', 1, ' ')
      end if

      ! Per bin stats ..
      if (do_drydep) then
        do ibin = 1, NBIN
          call addfld(trim(btndname(igroup, ibin))//'VD', 'm/s     ', 1,    'A', &
               trim(sname) // ' dry deposition velocity', phys_decomp)
        end do
      end if

    end do

    do igas = 1, NGAS
      icnst = icnst4gas(igas)

      ! Indicate that this is a constituent whose tendency could be changed by
      ! CARMA.
      lq_carma(icnst) = .true.
      gtndname(igas) = trim(cnst_name(icnst)) // 'TC'

      call addfld(gtndname(igas),     'kg/kg/s ', pver, 'A', &
           trim(cnst_name(icnst)) // ' CARMA tendency', phys_decomp)

      call addfld(trim(cnst_name(icnst))//'SI', 'ratio   ', pver, 'A', &
           trim(cnst_name(icnst)) // ' saturation wrt ice', phys_decomp)
      call addfld(trim(cnst_name(icnst))//'SL', 'ratio   ', pver, 'A', &
           trim(cnst_name(icnst)) // ' saturation wrt liquid', phys_decomp)
      call addfld(trim(cnst_name(icnst))//'EI', 'mol/mol ', pver, 'A', &
           trim(cnst_name(icnst)) // ' equilibrium vmr wrt ice', phys_decomp)
      call addfld(trim(cnst_name(icnst))//'EL', 'mol/mol ', pver, 'A', &
           trim(cnst_name(icnst)) // ' equilibrium vmr wrt liquid', phys_decomp)
      call addfld(trim(cnst_name(icnst))//'WT', '%       ', pver, 'A', &
           trim(cnst_name(icnst)) // ' weight percent aerosol composition', phys_decomp)
      
      call add_default(trim(cnst_name(icnst))//'SI', 1, ' ')
      call add_default(trim(cnst_name(icnst))//'SL', 1, ' ')
    end do
    
    if (carma_do_thermo) then
       call addfld('CRTT',     'K/s ', pver, 'A', ' CARMA temperature tendency', phys_decomp)
    end if
 
    ! Add fields for diagnostic fields, and make them defaults on the first tape.
    if (carma_do_substep) then
      call addfld('CRNSTEP',  '        ', pver, 'A', 'number of carma substeps', phys_decomp)
      call addfld('CRLNSTEP', '        ', pver, 'A', 'ln(number of carma substeps)', phys_decomp)
     
      call add_default('CRNSTEP', 1, ' ')
      call add_default('CRLNSTEP', 1, ' ')
    end if
    
    
    ! Set up the reference atmosphere that can be used for fixed initialization. This is
    ! an approximate atmospheric used to define average fall velocities, coagulation
    ! kernels, and growth parameters.
    if (carma_do_fixedinit) then

      ! NOTE: Reading the initial condtion file using the supplied routines must
      ! be done outside of masterproc, so does this in all threads before deciding
      ! if it will be used. The initial condition file is only opened on an initial run.
      if (is_first_step()) then      
        call carma_getT(carma_t_ref)
        if (carma%f_igash2o /= 0)    call carma_getH2O(carma_h2o_ref)
        if (carma%f_igash2So4 /= 0)  call carma_getH2SO4(carma_h2so4_ref)
      end if

      if (masterproc) then
        call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
        if (rc < 0) call endrun('carma_init::CARMA_Get failed.')
      
        if (do_print) write(LUNOPRT,*) ""
        if (do_print) write(LUNOPRT,*) "CARMA initializing to fixed reference state."
        if (do_print) write(LUNOPRT,*) ""
      
        ! For temperature, get the average temperature from reference temperature file
        ! if it exists or from the initial condition file if the reference temperature file
        ! doesn't exist.
        !
        ! NOTE: The reference temperature file will only be created for an inital run. It
        ! must already exist for a restart run.
        
        ! Does reference temperature file already exist?
        call getfil(carma_reftfile, locfn, iflag=1)
  
        inquire(file=locfn, exist=lexist)
  
        ! Read the reference temperature from the file.
        if (lexist) then
        
          ! Open the netcdf file.
          call wrap_open(trim(locfn), NF90_NOWRITE, ncid)
  
          ! Inquire about dimensions
          call wrap_inq_dimid(ncid, 'lev', dimid_lev)
          call wrap_inq_dimlen(ncid, dimid_lev, nlev)
          
          ! Does the number of levels match?
          if (nlev /= pver) then
            call endrun("carma_init::ERROR - Incompatible number of levels &
                 &in the CARMA reference temperature file ... " // trim(locfn))
          end if
  
          ! Get variable ID for reference temperature
          call wrap_inq_varid(ncid, 'T', vid_T)
  
          ! Read in the temperature data.
          call wrap_get_var_realx(ncid, vid_T, carma_T_ref)

          if (carma%f_igash2o /= 0) then
            ! Get variable ID for reference temperature
            call wrap_inq_varid(ncid, 'Q', vid_T)
  
            ! Read in the temperature data.
            call wrap_get_var_realx(ncid, vid_T, carma_h2o_ref)
          end if

          if (carma%f_igash2so4 /= 0) then
            ! Get variable ID for reference temperature
            call wrap_inq_varid(ncid, 'H2SO4', vid_T)
  
            ! Read in the temperature data.
            call wrap_get_var_realx(ncid, vid_T, carma_h2so4_ref)
          end if
          
          ! Close the file
          call wrap_close(ncid)
  
        ! Is this an initial or restart run?
        else if (is_first_step()) then

          if (do_print) write(LUNOPRT,*) ""
          if (do_print) write(LUNOPRT,*) 'Creating CARMA reference temperature file ... ', trim(locfn)
    
          ! Save the average into a file to be used for restarts.
          call CARMA_CreateRefTFile(carma, locfn, pref_mid(:) / 100._r8, &
               carma_t_ref(:), rc, refh2o=carma_h2o_ref(:), refh2so4=carma_h2so4_ref(:))
        else
          
          ! The file must already exist for a restart run.
          call endrun("carma_init::ERROR - Can't find the CARMA reference temperature file ... " // trim(carma_reftfile))

        end if

        ! Write out the values that are being used.
        if (do_print) write(LUNOPRT,*) ""
        if (do_print) write(LUNOPRT,1) "Level","Int P (Pa)","Mid P (Pa)","Mid T (K)"
        
        do iz = 1, pver
          if (do_print) write(LUNOPRT,2) iz, pref_edge(iz), pref_mid(iz), carma_t_ref(iz)
        end do
        if (do_print) write(LUNOPRT,2) iz, pref_edge(iz), 0.0_r8, 0.0_r8
        if (do_print) write(LUNOPRT,*) ""
      end if
      
#ifdef SPMD

      ! Communicate the settings to the other MPI tasks.
       call mpi_bcast(carma_t_ref,    pver,  MPI_REAL8, 0, mpicom, ier)
#endif
    end if


    ! Do a model specific initialization.
    call CARMA_InitializeModel(carma, lq_carma, rc)
    if (rc < 0) call endrun('carma_init::CARMA_InitializeModel failed.')
    
    return
  end subroutine carma_init


  !! Finalize (cleanup allocations) in the CARMA model.
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author  Chuck Bardeen
  !! @version October 2009
  subroutine carma_final
    implicit none
    
    integer           :: rc           ! CARMA return code
    integer           :: LUNOPRT      ! logical unit number for output
    logical           :: do_print     ! do print output?
            
    2 format(' carma_final: overall substepping statistics',/,&
           '    max nsubstep=',1F9.0,/,'    avg nsubstep=',1F9.2,/,&
           '    max nretry=',1F9.0,/,'    avg nretry=',1F10.4)

    ! Initialize the return code.
    rc = 0
    
    ! Output the end of run statistics for CARMA
    if (carma_do_substep) then
      if (masterproc) then
        call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
        if (rc < 0) call endrun('carma_final::CARMA_Get failed.')

        if (glob_nstep > 0) then
          if (do_print) write(LUNOPRT,2) glob_max_nsubstep, &
                                                     glob_nsubstep / glob_nstep, &
                                                     glob_max_nretry, &
                                                     glob_nretry / glob_nstep
        else
          if (do_print) write(LUNOPRT,2) glob_max_nsubstep, &
                                                     0., &
                                                     glob_max_nretry, &
                                                     0.
        end if
      end if
    end if
    
    
    ! Do a model specific initialization.
    call CARMA_Destroy(carma, rc)
    if (rc < 0) call endrun('carma_final::CARMA_Destroy failed.')
    
    return
  end subroutine carma_final


  !! Initialization that needs to be done prior to each timestep.
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine carma_timestep_init
    implicit none

    if (.not. carma_flag) return

    ! Reset the stats, so that they are per timestep values.
    step_max_nsubstep = 0._f
    step_max_nretry   = 0._f
    step_nstep        = 0._f
    step_nsubstep     = 0._f
    step_nretry       = 0._f
    
    return
  end subroutine carma_timestep_init


  !! Calculates the tendencies for all of the constituents handled by CARMA.
  !! To do this:
  !!
  !!  - a CARMASTATE object is created
  !!  - it is set to the current CAM state
  !!  - a new state is determined by CARMA
  !!  - the difference between these states is used to determine the tendencies
  !!  - statistics arecollected and reported
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! NOTE: Need to add code for getting/putting last fields into the physics
  !! buffer from substeping.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine carma_timestep_tend(state, cam_in, cam_out, ptend, dt, pbuf, dlf, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, ustar, obklen)
    use time_manager,     only: get_nstep, get_step_size, is_first_step
    use camsrfexch,       only: cam_in_t, cam_out_t
    use scamMod,          only: single_column
    use planck,           only: planckIntensity
    
    implicit none

    type(physics_state), intent(in)    :: state                 !! physics state variables
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    type(cam_out_t), intent(inout)     :: cam_out               !! cam output to surface models
    type(physics_ptend), intent(out)   :: ptend                 !! constituent tendencies
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    real(r8), intent(in), optional     :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
    real(r8), intent(inout), optional  :: rliq(pcols)           !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(out), optional    :: prec_str(pcols)       !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(out), optional    :: snow_str(pcols)       !! [Total] sfc flux of snow from stratiform   (m/s)
    real(r8), intent(out), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(out), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(in), optional     :: ustar(pcols)          !! friction velocity (m/s)
    real(r8), intent(in), optional     :: obklen(pcols)         !! Obukhov length [ m ]

    ! Local variables
    type(physics_state)   :: state_loc                              ! local physics state using wet mmr
    type(carma_type), pointer :: carma_ptr                          ! the carma state object
    type(carmastate_type) :: cstate                                 ! the carma state object
    integer               :: igroup                                 ! group index
    integer               :: ielem                                  ! element index
    integer               :: ielem_nd                               ! index of numder density element in group
    integer               :: ibin                                   ! bin index
    integer               :: igas                                   ! gas index
    integer               :: icol                                   ! column index
    integer               :: icnst                                  ! constituent index
    integer               :: icnst_q                                ! H2O constituent index
    integer               :: ncol                                   ! number of columns
    integer               :: rc                                     ! CARMA return code
    integer               :: cnsttype                               ! constituent type
    integer               :: maxbin                                 ! last prognostic bin
    real(r8)              :: spdiags(pcols, pver, NSPDIAGS)         ! CARMA step diagnostic output       
    real(r8)              :: gsdiags(pcols, pver, NGAS,   NGSDIAGS) ! CARMA gas diagnostic output       
    real(r8)              :: gpdiags(pcols, pver, NGROUP, NGPDIAGS) ! CARMA group diagnostic output 
    real(r8)              :: sbdiags(pcols, NBIN, NELEM,  NSBDIAGS) ! CARMA surface bin diagnostic output 
    real(r8)              :: bndiags(pcols, pver, NBIN, NELEM, NBNDIAGS) ! CARMA bin diagnostic output 
    real(r8)              :: newstate(pver)                         ! next state for a physics state field
    real(r8)              :: xc(pver)                               ! x center
    real(r8)              :: dx(pver)                               ! x width
    real(r8)              :: yc(pver)                               ! y center
    real(r8)              :: dy(pver)                               ! y width
    real(r8)              :: dz(pver)                               ! z width
    real(r8)              :: satice(pver)                           ! saturation wrt ice
    real(r8)              :: satliq(pver)                           ! saturation wrt liquid
    real(r8)              :: eqice(pver)                            ! equil vp wrt ice
    real(r8)              :: eqliq(pver)                            ! equil vp wrt liquid
    real(r8)              :: wtpct(pver)                            ! weight percent aerosol composition
    real(r8)              :: time                                   ! the total elapsed time (s)
    real(r8)              :: dlat                                   ! latitude spacing
    real(r8)              :: r(NBIN)                                ! particle radius (cm)
    real(r8)              :: rmass(NBIN)                            ! particle mass (g)
    real(r8)              :: rrat(NBIN)                             ! particle maximum radius ratio ()
    real(r8)              :: arat(NBIN)                             ! particle area ration ()
    real(r8)              :: rhoelem                                ! element density (g)
    real(r8)              :: nd(pver)                               ! number density (cm-3)
    real(r8)              :: ad(pver)                               ! area density (um2/cm3)
    real(r8)              :: md(pver)                               ! mass density (g cm-3)
    real(r8)              :: mr(pver)                               ! mass mixing ratio (kg/kg)
    real(r8)              :: re(pver)                               ! effective radius (um)
    real(r8)              :: rm(pver)                               ! Mitchell effective radius (um)
    real(r8)              :: ex(pver)                               ! extinction (km-1)
    real(r8)              :: od(pver)                               ! optical depth
    real(r8)              :: re2(pver)                              ! N(r)*r^2 (cm2)
    real(r8)              :: re3(pver)                              ! N(r)*r^3 (cm3)
    real(r8)              :: pa(pver)                               ! Projected Area (cm2)
    real(r8)              :: ar(pver)                               ! Area Ratio 
    real(r8)              :: vm(pver)                               ! Massweighted fall velocity (cm2)
    real(r8)              :: jn(pver)                               ! nucleation (cm-3)
    real(r8)              :: numberDensity(pver)                    ! number density (cm-3)
    real(r8)              :: nucleationRate(pver)                   ! nucleation rate (cm-3 s-1)
    real(r8)              :: extinctionCoefficient(pver)            ! extinction coefficient (cm2)
    real(r8)              :: dd                                     ! dry deposition (kg/m2)
    real(r8)              :: vd                                     ! dry deposition velocity (cm/s)
    real(r8)              :: vf(pverp)                              ! fall velocity (cm/s)
    real(r8)              :: dtpart(pver)                           ! delta particle temperature (K)
    real(r8), pointer, dimension(:, :) :: t_ptr                     ! last temperature
    real(r8), pointer, dimension(:, :) :: gc_ptr                    ! last gas mmr
    real(r8), pointer, dimension(:, :) :: sati_ptr                  ! last saturation wrt ice
    real(r8), pointer, dimension(:, :) :: satl_ptr                  ! last saturation wrt liquid
    real(r8), pointer, dimension(:, :, :) :: su_ptr                 ! shortwave flux up (W/m2)
    real(r8), pointer, dimension(:, :, :) :: sd_ptr                 ! shortwave flux down (W/m2)
    real(r8), pointer, dimension(:, :, :) :: lu_ptr                 ! longwave flux up (W/m2)
    real(r8), pointer, dimension(:, :, :) :: ld_ptr                 ! longwave flux down (W/m2)
    real(r8), pointer, dimension(:,:) :: tnd_qsnow    ! external tendency on snow mass (kg/kg/s)
    real(r8), pointer, dimension(:,:) :: tnd_nsnow    ! external tendency on snow number(#/kg/s)
    real(r8), pointer, dimension(:,:) :: re_ice       ! ice effective radius (m)
    integer               :: lchnk                                  ! chunk identifier
    real(r8)              :: coremmr(pver)
    real(r8)              :: ttlmmr(pver)
    integer               :: iz
    real(r8)              :: cldfrc(pver)                           ! cloud fraction [fraction]
    real(r8)              :: rhcrit(pver)                           ! relative humidity for onset of liquid clouds [fraction]
    real(r8)              :: lndram                                 ! land aerodynamical resistance [s/m]
    real(r8)              :: lndfv                                  ! surface friction velocity from land [m/s]
    real(r8)              :: ocnram                                 ! ocean aerodynamical resistance [s/m]
    real(r8)              :: ocnfv                                  ! surface friction velocity from ocean [m/s]
    real(r8)              :: iceram                                 ! sea ice aerodynamical resistance [s/m]
    real(r8)              :: icefv                                  ! surface friction velocity from sea ice [m/s]
    real(r8)              :: radint(pver,NWAVE)                     ! radiative intensity (W/m2/sr/cm)
    real(kind=f)          :: dwave(NWAVE)                           ! the wavelengths widths (cm)
    real(kind=f)          :: wave(NWAVE)                            ! the center wavelengths (cm)

    integer               :: max_nsubstep
    real(kind=f)          :: max_nretry
    real(kind=f)          :: nstep
    integer               :: nsubstep
    real(kind=f)          :: nretry
    real(kind=f)          :: zsubsteps(pver)
    logical               :: is_cloud                               ! is the group a cloud?
    logical               :: is_ice                                 ! is the group ice?
    integer               :: ienconc
    logical               :: grp_do_drydep                          ! is dry depostion enabled for group?
    logical               :: do_drydep                              ! is dry depostion enabled?
    logical               :: do_fixedinit                           ! do initialization from reference atm?
    logical               :: do_detrain                             ! do convective detrainment?
    integer               :: iwvl
    real(r8), parameter   :: zzocen = 0.0001_r8                     ! Ocean aerodynamic roughness length [m]
    real(r8), parameter   :: zzsice = 0.0400_r8                     ! Sea ice aerodynamic roughness length [m]

    
    ! Initialize the return code.
    rc = 0

    ! Initialize the output tendency structure.
    call physics_ptend_init(ptend,state%psetcols,'CARMA', ls=carma_do_thermo, lq=lq_carma)
    
    if (present(prec_sed)) prec_sed(:) = 0._f
    if (present(snow_sed)) snow_sed(:) = 0._f
    if (present(prec_str)) prec_str(:) = 0._f
    if (present(snow_str)) snow_str(:) = 0._f
    
    if (.not. carma_flag) return

    ! Determine the current time in seconds.
    time = dt * get_nstep() - 1
    
    ! The CARMA interface assumes that mass mixing ratios are relative to a
    ! wet atmosphere, so convert any dry mass mixing ratios to wet.
    call physics_state_copy(state, state_loc)
    call set_dry_to_wet(state_loc)
    
    spdiags(:, :, :)       = 0.0_r8
    gpdiags(:, :, :, :)    = 0.0_r8
    gsdiags(:, :, :, :)    = 0.0_r8
    sbdiags(:, :, :, :)    = 0.0_r8
    bndiags(:, :, :, :, :) = 0.0_r8
    
    ! Find the constituent index for water vapor.
    call cnst_get_ind('Q', icnst_q)
    
    ! Get pointers into pbuf ...
    lchnk = state_loc%lchnk
    
    call pbuf_get_field(pbuf, ipbuf4t, t_ptr)
    
    ! If doing particle heating, then get pointers to the spectral flux data provided
    ! by the radiation code in the physics buffer.
    !
    ! NOTE: RRTMG can now be done in a subset of all levels, rather than the full
    ! model height. Any implications for this code have not been worked out.
    if (carma_do_pheat) then
      call pbuf_get_field(pbuf, pbuf_get_index("SU"), su_ptr)
      call pbuf_get_field(pbuf, pbuf_get_index("SD"), sd_ptr)
      call pbuf_get_field(pbuf, pbuf_get_index("LU"), lu_ptr)
      call pbuf_get_field(pbuf, pbuf_get_index("LD"), ld_ptr)
    end if
    
    ! Cloud ice pbuf fields
    if (carma_do_cldice) then
      call pbuf_get_field(pbuf, pbuf_get_index("TND_QSNOW"), tnd_qsnow)
      call pbuf_get_field(pbuf, pbuf_get_index("TND_NSNOW"), tnd_nsnow)
      call pbuf_get_field(pbuf, pbuf_get_index("RE_ICE"), re_ice)
    end if


    ! Create a CARMASTATE object which contains state information about one
    ! column of the atmosphere.
    carma_ptr => carma


    ! If initializing CARMASTATE from a reference state, do it before entering the main
    ! loop.
    !
    call CARMA_Get(carma, rc, do_fixedinit=do_fixedinit, do_drydep=do_drydep)
    if (rc < 0) call endrun('carma_timestep_tend::CARMA_Get failed.')
    
    if (do_fixedinit) then
    
      ! The latitude and longitude are arbitrary, but the dimensions need to be correct.
      xc = 255._r8
      yc = 40._r8

      ! Assume resolution is 64x128.
      if (single_column) then
        dx = 360._r8 / 128._r8
        dy = 180._r8 / 64._r8
      else
      
        ! Calculate the x and y coordinates, in degrees latitude and longitude.
        dx = 360._r8 / plon
        dy = 180._r8 / (plat-1)
      end if

      call CARMASTATE_CreateFromReference(cstate, &
                         carma_ptr, &
                         time, &
                         dt, &
                         pver, &
                         I_HYBRID, &
                         I_LL, &
                         40._r8, &
                         255._r8, &
                         xc, &
                         dx, &
                         yc, &
                         dy, &
                         pref_mid_norm, &
                         pref_edge/psurf_ref, &
                         pref_mid(:), &
                         pref_edge(:), &
                         carma_t_ref(:), &
                         rc, &
                         qh2o=carma_h2o_ref, &
                         qh2so4=carma_h2so4_ref)
      if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_CreateFromReference failed.')
    end if


    ! Process each column.
    do icol = 1, state_loc%ncol
    
      ! Haven't figured out how to get dimensions for single column. Perhaps should change
      ! CARMA to work with area rather than dx and dy. For now, just hack something.
      xc(:) = state_loc%lon(icol) / DEG2RAD
      yc(:) = state_loc%lat(icol) /  DEG2RAD

      ! Assume resolution is 64x128.
      if (single_column) then
        dx = 360._r8 / 128._r8
        dy = 180._r8 / 64._r8
      else
      
        ! Caclulate the x and y coordinates, in degrees latitude and longitude.
        dx(:) = 360._r8 / plon

        dlat = 180._r8 / (plat-1)

        ! The pole points need special treatment, since the point is not the
        ! center of the grid box.
        !
        ! In single column mode there is just one latitude, so make it global.
        if (abs(state_loc%lat(icol) /  DEG2RAD) >= (90._r8 - (90._r8 / (plat-1)))) then

          ! Nudge yc toward the equator.
          yc(:) = yc(:) - sign(0.25_r8,state_loc%lat(icol)) * dlat
          
          dy(:) = dlat / 2._r8
        else
          dy(:) = dlat
        endif
      end if

      if (is_first_step()) then
        t_ptr(icol,:) = state_loc%t(icol,:)
      end if
      
      ! For particle heating, need to get the incoming radiative intensity from
      ! the radiation code.
      !
      ! The radiation code can optionally provide the flux up and down per band in W/m2,
      ! when the compute_spectral_flux namelist variable is provided to the radiation. This
      ! data needs to be scaled to a radiative intensity by assuming it is isotrotropic.
      radint(:,:) = 0._f
      
      if (carma_do_pheat) then
        call CARMA_Get(carma, rc, dwave=dwave, wave=wave)
        if (rc < 0) call endrun('carma_timestep_tend::CARMA_Get failed.')
        
        ! CARMA may run before the radiation code for the very first time step.
        ! In that case, the lu, ld, su and sd values are NaN. NaN will crash
        ! the model, so instead substitute an approximation that is roughly a
        ! nighttime (su=sd=0) with a black body temperature of the grid point
        ! temperature (lu=ld=B(T)).
        !
        ! NOTE: planckIntensity is in erg/cm2/s/sr/cm and lu is in W/m2,
        ! so some conversion factors are needed.
        if (is_first_step()) then
          su_ptr(icol, :, :) = 0._r8
          sd_ptr(icol, :, :) = 0._r8

          do iwvl = 1, nlwbands
            do iz = 1, pver
              lu_ptr(icol, iz, iwvl) = planckIntensity(wave(iwvl), state_loc%t(icol, iz)) / 1e7_f * 1e4_f * dwave(iwvl) * PI
            end do
            lu_ptr(icol, pverp, iwvl) = lu_ptr(icol, pver, iwvl)
            
            ld_ptr(icol, 2:pverp, iwvl) = lu_ptr(icol, 1:pver, iwvl)
            ld_ptr(icol, 1, iwvl) = lu_ptr(icol, 2, iwvl)
          end do
        end if

        
        do iwvl = 1, nlwbands
          radint(:, iwvl) = (lu_ptr(icol, 2:, iwvl) + ld_ptr(icol, :pver, iwvl)) / 2._r8 / PI / dwave(iwvl)
        end do
        
        do iwvl = 1, nswbands
          radint(:, nlwbands+iwvl) = (su_ptr(icol, 2:, iwvl) + sd_ptr(icol, :pver, iwvl)) / 2._r8 / PI / dwave(nlwbands+iwvl)
        end do
      end if

      call CARMASTATE_Create(cstate, &
                             carma_ptr, &
                             time, &
                             dt, &
                             pver, &
                             I_HYBRID, &
                             I_LL, &
                             state_loc%lat(icol) / DEG2RAD, &
                             state_loc%lon(icol) / DEG2RAD, &
                             xc, &
                             dx, &
                             yc, &
                             dy, &
                             pref_mid_norm, &
                             pref_edge/psurf_ref, &
                             state_loc%pmid(icol, :), &
                             state_loc%pint(icol, :), &
                             state_loc%t(icol, :), &
                             rc, &
                             qh2o=state_loc%q(icol, :, icnst_q), &
                             told=t_ptr(icol, :), &
                             radint=radint)
      if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_Create failed.')


      ! Store information about the CARMA particles.

      ! For prognostic groups, the mass of the particles for each bin is stored as
      ! a unique constituent within CAM.
      do ielem = 1, NELEM
        call CARMAELEMENT_Get(carma,  ielem, rc, igroup=igroup)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAELEMENT_Get failed.')
        
        call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, maxbin=maxbin)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAGROUP_Get failed.')
        
        if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then

          ! For prognostic groups, set the bin from the corresponding constituent.
          do ibin = 1, NBIN

            ! Bins past maxbin are treated as diagnostic even if the group
            ! is prognostic and thus are not advected in the parent model.
            if (ibin <= maxbin) then
              call CARMASTATE_SetBin(cstate, ielem, ibin, state_loc%q(icol, :, icnst4elem(ielem, ibin)), rc)
              if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_SetBin failed.')
            else
              newstate(:) = 0._f
              
              call CARMASTATE_SetBin(cstate, ielem, ibin, newstate, rc)
              if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_SetBin failed.')
            end if
          end do
        end if
      end do

      ! Store information about CARMA gases.
      do igas = 1, NGAS
        call pbuf_get_field(pbuf, ipbuf4gas(igas), gc_ptr)
        call pbuf_get_field(pbuf, ipbuf4sati(igas), sati_ptr)
        call pbuf_get_field(pbuf, ipbuf4satl(igas), satl_ptr)

        ! Handle the initial case where we don't have last values.
        if (is_first_step()) then
          gc_ptr(icol,:) = state_loc%q(icol, :, icnst4gas(igas))
          sati_ptr(icol,:) = -1._f
          satl_ptr(icol,:) = -1._f
        end if

        call CARMASTATE_SetGas(cstate, igas, state_loc%q(icol, :, icnst4gas(igas)), rc, &
          mmr_old=gc_ptr(icol,:), satice_old=sati_ptr(icol,:), satliq_old=satl_ptr(icol,:))
        if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_SetGas failed.')
      end do


      call CARMA_DiagnoseBins(carma, cstate, state_loc, pbuf, icol, dt, rc, rliq=rliq, prec_str=prec_str, snow_str=snow_str)
      if (rc < 0) call endrun('carma_timestep_tend::CARMA_DiagnoseBins failed.') 
      
      
      ! If the model supports detraining of condensed water from convection, then pass
      ! along the condensed H2O.
      call CARMA_Get(carma, rc, do_detrain=do_detrain)
      if (rc < 0) call endrun('CARMA_Detrain::CARMA_Get failed.')

      if (do_detrain) then
        call CARMA_Detrain(carma, cstate, cam_in, dlf, state_loc, icol, dt, rc, rliq=rliq, prec_str=prec_str, &
              snow_str=snow_str, tnd_qsnow=tnd_qsnow, tnd_nsnow=tnd_nsnow)
        if (rc < 0) call endrun('carma_timestep_tend::CARMA_Detrain failed.')
      end if
      

      ! Now that detrainment has happened, determine the cloud fractions.
      ! These will be used to scale the cloud amount to go from gridbox average to in-cloud
      ! values and back.
      !
      ! For the cirrus model, assume the cloud fraction is just the ice cloud fraction.
      call CARMA_CloudFraction(carma, cstate, cam_in, state_loc, icol, cldfrc, rhcrit, rc)
      if (rc < 0) call endrun('carma_timestep_tend::carma_CloudFraction failed.')

      ! A fixed value for rhcrit can be specified in the namelist rather than using the
      ! one from the cloud fraction.
      if (carma_rhcrit /= 0._f) then
        rhcrit(:) = carma_rhcrit
      end if
      
      
      ! For dry deposition, provide a surface friction velocity and an aerodynamic
      ! resistance for each of the land surface types. The values for the land come
      ! from the land model, but those for ocean and sea ice need to be calculated.
      if (do_drydep) then
      
        ! Land
        lndfv = cam_in%fv(icol)
        lndram = cam_in%ram1(icol)
        
        ! Ocean
        ocnfv  = ustar(icol)
        ocnram = 0._r8
        if (cam_in%ocnfrac(icol) > 0._r8) then
          call CARMA_calcram(ocnfv, &
                             zzocen, &
                             state_loc%pdel(icol, pver), &
                             state_loc%pmid(icol, pver), &
                             state_loc%t(icol, pver), &
                             obklen(icol), &
                             ocnram)
        end if

        ! Sea Ice
        icefv  = ustar(icol)
        iceram = 0._r8
        if (cam_in%icefrac(icol) > 0._r8) then
          call CARMA_calcram(ocnfv, &
                             zzocen, &
                             state_loc%pdel(icol, pver), &
                             state_loc%pmid(icol, pver), &
                             state_loc%t(icol, pver), &
                             obklen(icol), &
                             iceram)
        end if
      end if
    
    
      ! Advance the microphysics one timestep.
      call CARMASTATE_Step(cstate, rc, cldfrc=cldfrc, rhcrit=rhcrit, &
           lndfv=lndfv, ocnfv=ocnfv, icefv=icefv, lndram=lndram, &
           ocnram=ocnram, iceram=iceram, lndfrac=cam_in%landfrac(icol), &
           ocnfrac=cam_in%ocnfrac(icol), icefrac=cam_in%icefrac(icol))
      if (rc < 0) call endrun('carma_timestep_tend::CARMA_Step failed.') 
            

      ! Get the results for the CARMA particles.

      ! For diagnostic groups, a special routine needs to be called to determine how
      ! bins affect the bulk state, since there is not an individual constituent for
      ! each bin.
      !
      ! NOTE: To work around an XL Fortran compiler bug, the optional arguments can only
      ! be passed when defined.
      if (present(rliq)) then
        call CARMA_DiagnoseBulk(carma, cstate, cam_out, state_loc, pbuf, ptend, icol, dt, rc, &
          rliq=rliq, prec_str=prec_str, snow_str=snow_str, prec_sed=prec_sed, &
          snow_sed=snow_sed, tnd_qsnow=tnd_qsnow, tnd_nsnow=tnd_nsnow, re_ice=re_ice)
      else
        call CARMA_DiagnoseBulk(carma, cstate, cam_out, state_loc, pbuf, ptend, icol, dt, rc)
      end if
      if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_DiagnoseBulk failed.')


      ! Calculate the group statistics for all elements.
      dz(:) = state_loc%zi(icol, 1:pver) - state_loc%zi(icol, 2:pverp)
            
      do ielem = 1, NELEM
      
        call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAELEMENT_Get failed.')
        
        call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, r=r, rmass=rmass, maxbin=maxbin, &
               is_cloud=is_cloud, is_ice=is_ice, do_drydep=grp_do_drydep, rrat=rrat, arat=arat)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAGROUP_Get failed.')
      
        ! Intialize the group totals
        nd(:)  = 0.0_r8
        ad(:)  = 0.0_r8
        md(:)  = 0.0_r8
        mr(:)  = 0.0_r8
        re(:)  = 0.0_r8
        rm(:)  = 0.0_r8
        ex(:)  = 0.0_r8
        od(:)  = 0.0_r8
        re2(:) = 0.0_r8
        re3(:) = 0.0_r8
        jn(:)  = 0.0_r8
        pa(:)  = 0.0_r8
        vm(:)  = 0.0_r8
      
        do ibin = 1, NBIN
          call CARMASTATE_GetBin(cstate, ielem, ibin, newstate(:), rc, &
                 numberDensity=numberDensity, nucleationRate=nucleationRate, surface=dd, vd=vd, vf=vf, dtpart=dtpart)
          if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_GetBin failed.')
          
          ! For prognostic groups, set the tendency from the corresponding constituents.
          if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then
          
            ! Bins past maxbin are treated as diagnostic even if the group
            ! is prognostic and thus are not advected in the paerent model.
            if (ibin <= maxbin) then
    
              icnst = icnst4elem(ielem, ibin)
            
              ! Update the consituent tendency.
              ptend%q(icol, :, icnst) = (newstate(:) - state_loc%q(icol, :, icnst)) / dt
              
              if (grp_do_drydep) then
                sbdiags(icol, ibin, ielem, SBDIAGS_DD) = dd / dt
                sbdiags(icol, ibin, ielem, SBDIAGS_VD) = - vd / 100._r8
              end if
            end if
          end if
          
          ! Calculate the total densities.
          !
          ! NOTE: Convert AD to um2/cm3.
          if (numberDensity(1) /= CAM_FILL) then
            nd(:)  = nd(:)  + numberDensity(:)
            re2(:) = re2(:) + numberDensity(:) * ((r(ibin)*rrat(ibin))**2)
            re3(:) = re3(:) + numberDensity(:) * ((r(ibin)*rrat(ibin))**3)
            ad(:)  = ad(:)  + numberDensity(:) * 4.0_r8 * PI * (r(ibin)**2) * 1.0e8_r8
            md(:)  = md(:)  + numberDensity(:) * rmass(ibin)
            mr(:)  = mr(:)  + newstate(:)
            pa(:)  = pa(:)  + numberDensity(:) * PI * ((r(ibin) * rrat(ibin))**2) * arat(ibin)
            vm(:)  = vm(:)  + numberDensity(:) * rmass(ibin) * vf(2:) / 100._f
        
            ! Calculate the optical depth and extinction.
            !
            ! NOTE: Assume Qext = 2 for optical depth. This can be pulled out of CARMA
            ! mie claculations later.
            !
            ! Convert extinction coefficient to km-1.
            extinctionCoefficient(:) = 2.0_r8 * PI * (r(ibin)**2)
            ex(:) = ex(:) + numberDensity(:) * extinctionCoefficient(:) * 1e5_r8
            od(:) = od(:) + numberDensity(:) * extinctionCoefficient(:) * dz(:) * 100._r8
          end if

          ! Particle temperatures from particle heating.
          if (carma_do_pheat) then
            bndiags(icol, :, ibin, ielem, BNDIAGS_TP) = dtpart(:)
          end if

          if (nucleationRate(1) /= CAM_FILL) then
            jn(:)  = jn(:)  + nucleationRate(:)
          end if
        end do
      
        ! If this is the number element for the group, then write out the
        ! statistics.
        if (numberDensity(1) /= CAM_FILL) then
          
          ! Calculate the effective radius (total volume / total area). Places
          ! with no surface area will cause NaN values.
          !
          ! NOTE: Convert RE to um.
          where (re2(:) > 0.0_r8)
            re(:) = (re3(:) / re2(:)) * 1e4_r8
            rm(:) = (3._r8 / 4._r8) * (md(:)  / (0.917_r8 * pa(:))) * 1e4_r8
            ar(:) = pa(:) / PI / re2(:)
          end where

          where (md(:) > 0.0_r8)
            vm(:) = vm(:) / md(:)
          end where

          ! Store the statistics.
          
          ! Gridbox average
          gpdiags(icol, :, igroup, GPDIAGS_ND) = nd 
          gpdiags(icol, :, igroup, GPDIAGS_AD) = ad
          gpdiags(icol, :, igroup, GPDIAGS_MD) = md
          gpdiags(icol, :, igroup, GPDIAGS_RE) = re
          gpdiags(icol, :, igroup, GPDIAGS_RM) = rm
          gpdiags(icol, :, igroup, GPDIAGS_MR) = mr
          gpdiags(icol, :, igroup, GPDIAGS_EX) = ex
          gpdiags(icol, :, igroup, GPDIAGS_OD) = od
          gpdiags(icol, :, igroup, GPDIAGS_VM) = vm
          gpdiags(icol, :, igroup, GPDIAGS_PA) = pa
          gpdiags(icol, :, igroup, GPDIAGS_AR) = ar
                  
          if (nucleationRate(1) /= CAM_FILL) then
            gpdiags(icol, :, igroup, GPDIAGS_JN) = jn
          end if
        end if
      end do

      
      ! Get the results for the CARMA gases.
      do igas = 1, NGAS
        call pbuf_get_field(pbuf, ipbuf4gas(igas), gc_ptr)
        call pbuf_get_field(pbuf, ipbuf4sati(igas), sati_ptr)
        call pbuf_get_field(pbuf, ipbuf4satl(igas), satl_ptr)

        call CARMASTATE_GetGas(cstate, igas, newstate(:), rc, satice=satice, satliq=satliq, &
         eqice=eqice, eqliq=eqliq, wtpct=wtpct)
        if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_GetGas failed.')
          
        icnst = icnst4gas(igas)

        ptend%q(icol, :, icnst) = (newstate(:) - state_loc%q(icol, :, icnst)) / dt
 
        gsdiags(icol, :, igas, GSDIAGS_SI) = satice(:) 
        gsdiags(icol, :, igas, GSDIAGS_SL) = satliq(:)
        gsdiags(icol, :, igas, GSDIAGS_EI) = eqice(:) 
        gsdiags(icol, :, igas, GSDIAGS_EL) = eqliq(:)
        gsdiags(icol, :, igas, GSDIAGS_WT) = wtpct(:)
        
        ! Store the values needed for substepping in the physics buffer.
        gc_ptr(icol,:)    = newstate(:)
        sati_ptr(icol, :) = satice(:)
        satl_ptr(icol, :) = satliq(:)
      end do

     
      ! Get the results for temperature.
      call CARMASTATE_GetState(cstate, rc, t=newstate(:))
      if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_GetState failed.')
          
      ! Store the values needed for substepping in the physics buffer.
      t_ptr(icol,:) = newstate(:)

      if (carma_do_thermo) then 
        ptend%s(icol, :) = (newstate(:) - state_loc%t(icol, :)) * cpair / dt
      endif
      
      
      ! Get the substepping statistics
      if (carma_do_substep) then
        call CARMASTATE_Get(cstate, rc, zsubsteps=zsubsteps)
        if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_Get failed.')
      
        spdiags(icol, :, SPDIAGS_NSTEP) = zsubsteps(:)
        spdiags(icol, :, SPDIAGS_LNSTEP) = log(zsubsteps(:))
      end if
    end do
    
    
    ! Report substep diagnostics
    if (carma_do_substep) then
      call CARMASTATE_Get(cstate, rc, max_nsubstep=max_nsubstep, max_nretry=max_nretry, &
        nstep=nstep, nsubstep=nsubstep, nretry=nretry)
      if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_Get failed.')

!$OMP CRITICAL
      step_max_nsubstep = max(step_max_nsubstep, real(max_nsubstep, f))
      step_max_nretry   = max(step_max_nretry, max_nretry)
    
      step_nstep        = step_nstep    + nstep
      step_nsubstep     = step_nsubstep + real(nsubstep, f)
      step_nretry       = step_nretry   + nretry
!$OMP END CRITICAL
    end if
     
    ! The CARMASTATE object is no longer needed.
    call CARMASTATE_Destroy(cstate, rc)
    if (rc < 0) call endrun('carma_timestep_tend::CARMASTATE_Destroy failed.')
    
     
    ! Output diagnostic fields.
    call carma_output_diagnostics(state_loc, ptend, gpdiags, sbdiags, gsdiags, spdiags, bndiags)

  end subroutine carma_timestep_tend
  
  
  subroutine carma_accumulate_stats()
    implicit none
    
    integer               :: istat
    integer               :: rc
    real(kind=f)          :: wrk
    integer               :: LUNOPRT              ! logical unit number for output
    logical               :: do_print             ! do print output?

    !  Define formats
    1 format(' carma: max nsubstep=',1F9.0,3x,'avg nsubstep=',1F9.2,3x,'max nretry=',1F9.0,3x,'avg nretry=',1F10.4)    

    if (carma_do_substep) then

      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun('carma_init::CARMA_Get failed.')

#ifdef SPMD
      call mpi_allreduce(step_max_nsubstep, wrk, 1, mpir8, mpi_max, mpicom, istat)
      if( istat /= MPI_SUCCESS ) then
         if (do_print) write(LUNOPRT,*) 'carma_timestep_tend: MPI_ALLREDUCE for max_nsubstep failed; error = ',istat
         call endrun
      end if
      step_max_nsubstep = wrk
      glob_max_nsubstep = max(glob_max_nsubstep, wrk)
      
      call mpi_allreduce(step_max_nretry, wrk, 1, mpir8, mpi_max, mpicom, istat)
      if( istat /= MPI_SUCCESS ) then
         if (do_print) write(LUNOPRT,*) 'carma_timestep_tend: MPI_ALLREDUCE for max_nsubstep failed; error = ',istat
         call endrun
      end if
      step_max_nretry = wrk
      glob_max_nretry = max(glob_max_nretry, wrk)
  
      call mpi_allreduce(step_nstep, wrk, 1, mpir8, mpi_sum, mpicom, istat)
      if( istat /= MPI_SUCCESS ) then
         if (do_print) write(LUNOPRT,*) 'carma_timestep_tend: MPI_ALLREDUCE for nstep failed; error = ',istat
         call endrun
      end if
      step_nstep = wrk
      glob_nstep = glob_nstep + wrk
  
      call mpi_allreduce(step_nsubstep, wrk, 1, mpir8, mpi_sum, mpicom, istat)
      if( istat /= MPI_SUCCESS ) then
         if (do_print) write(LUNOPRT,*) 'carma_timestep_tend: MPI_ALLREDUCE for nsubstep failed; error = ',istat
         call endrun
      end if
      step_nsubstep = wrk
      glob_nsubstep = glob_nsubstep + wrk
  
      call mpi_allreduce(step_nretry, wrk, 1, mpir8, mpi_sum, mpicom, istat)
      if( istat /= MPI_SUCCESS ) then
         if (do_print) write(LUNOPRT,*) 'carma_timestep_tend: MPI_ALLREDUCE for nretry failed; error = ',istat
         call endrun
      end if
      step_nretry = wrk
      glob_nretry = glob_nretry + wrk
#else

      ! For single CPU or OMP, just set the globals directly.
      glob_max_nsubstep = max(glob_max_nsubstep, step_max_nsubstep)
      glob_max_nretry   = max(glob_max_nretry, step_max_nretry)
      glob_nstep        = glob_nstep + step_nstep
      glob_nsubstep     = glob_nsubstep + step_nsubstep
      glob_nretry       = glob_nretry + step_nretry
      
#endif

      if (masterproc) then
        if (step_nstep > 0) then
          if (do_print) write(LUNOPRT,1) step_max_nsubstep, &
                                                     step_nsubstep / step_nstep, &
                                                     step_max_nretry, &
                                                     step_nretry / step_nstep
        else
          if (do_print) write(LUNOPRT,1) step_max_nsubstep, &
                                                     0., &
                                                     step_max_nretry, &
                                                     0.
        end if
      end if
    end if

  end subroutine carma_accumulate_stats


  !! Set initial mass mixing ratios of constituents, if nothing is specifed
  !! in the initial conditions file.
  !!
  !! NOTE: This call is part of the CAM Physics Interface
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine carma_init_cnst(name, q, gcid)
    implicit none

    character(len=*), intent(in) :: name               !! constituent name
    real(r8), intent(out)        :: q(:,:)             !! mass mixing ratio (gcol, lev)
    integer, intent(in)          :: gcid(:)            !! global column id
    
    integer                      :: igroup             ! group index
    integer                      :: ielem              ! element index
    integer                      :: ibin               ! bin index
    integer                      :: icnst              ! constituent index
    integer                      :: rc                 ! CARMA return code
    integer                      :: cnsttype           ! constituent type
    integer                      :: maxbin             ! last prognostic bin

    ! Initialize the return code.
    rc = 0
    
    ! Determine the element an bin for the particle
    do ielem = 1, NELEM
      do ibin = 1, NBIN
      
        call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAELEMENT_Get failed.')
        
        call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, maxbin=maxbin)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAGROUP_Get failed.')
        
        if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then

          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the paerent model.
          if (ibin <= maxbin) then
    
            icnst = icnst4elem(ielem, ibin)
            
            if (cnst_name(icnst) == name) then
    
              ! By default, initialize all constituents to 0.
              q(:, :) = 0.0_r8
              
              call CARMA_InitializeParticle(carma, ielem, ibin, q, gcid, rc)
              if (rc < 0) call endrun('carma_init_cnst::CARMA_InitializeParticle failed.')
            end if
          end if
        end if
      end do
    end do
    
    ! NOTE: There is currently no initialization for gases, but it could be
    ! added here.
    
    return
  end subroutine carma_init_cnst


  !! Outputs tracer tendencies and diagnositc fields to the history files.
  !! All the columns in the chunk should be output at the same time.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine carma_output_diagnostics(state, ptend, gpdiags, sbdiags, gsdiags, spdiags, bndiags)
    use cam_history, only: outfld

    implicit none

    type(physics_state), intent(in)   :: state          !! Physics state variables - before CARMA
    type(physics_ptend), intent(in)   :: ptend          !! indivdual parameterization tendencies
    real(r8), intent(in), dimension(pcols, pver, NGROUP, NGPDIAGS) :: gpdiags  !! CARMA group diagnostic output       
    real(r8), intent(in), dimension(pcols, NBIN, NELEM,  NSBDIAGS) :: sbdiags  !! CARMA surface bin diagnostic output 
    real(r8), intent(in), dimension(pcols, pver, NGAS,   NGSDIAGS) :: gsdiags  !! CARMA gas diagnostic output       
    real(r8), intent(in), dimension(pcols, pver, NSPDIAGS)         :: spdiags  !! CARMA step diagnostic output       
    real(r8), intent(in), dimension(pcols, pver, NBIN, NELEM, NBNDIAGS) :: bndiags !! CARMA bin diagnostic output 

    ! Local variables
    integer           :: igroup    ! group index
    integer           :: ielem     ! element index
    integer           :: ibin      ! bin index
    integer           :: igas      ! gas index
    integer           :: ienconc   ! element index for group's concentration element
    integer           :: icnst     ! constituent index
    integer           :: lchnk     ! chunk identifier
    integer           :: ncol      ! number of columns
    integer           :: rc        ! CARMA return code
    character(len=8)  :: sname     ! short (CAM) name
    integer           :: cnsttype  ! constituent type
    integer           :: maxbin    ! last prognostic bin
    logical           :: is_cloud  ! is the group a cloud?
    logical           :: do_drydep ! is dry deposition enabled?
       
    ! Initialize the return code.
    rc = 0
    
    ! Check each column int the chunk.
    lchnk = state%lchnk
    ncol  = state%ncol

    ! Output step diagnostics.
    if (carma_do_substep) then
      call outfld('CRNSTEP',  spdiags(:, :, SPDIAGS_NSTEP), pcols, lchnk) 
      call outfld('CRLNSTEP', spdiags(:, :, SPDIAGS_LNSTEP), pcols, lchnk) 
    end if

    ! Output the particle tendencies.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
      
        call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAELEMENT_Get failed.')
        
        call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, maxbin=maxbin, do_drydep=do_drydep)
        if (rc < 0) call endrun('carma_timestep_tend::CARMAGROUP_Get failed.')
        
        if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then

          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the paerent model.
          if (ibin <= maxbin) then
    
            icnst = icnst4elem(ielem, ibin)
  
            call outfld(trim(etndname(ielem, ibin))//'TC', ptend%q(:, :, icnst), pcols, lchnk) 

            if (do_drydep) then
              call outfld(trim(etndname(ielem, ibin))//'DD', sbdiags(:, ibin, ielem, SBDIAGS_DD), pcols, lchnk)
            end if

            if (carma_do_pheat) then
            
              ! Only specified for the number density element of the group.
              if (bndiags(1, 1, ibin, ielem, BNDIAGS_TP) /= CAM_FILL) then
                call outfld(trim(etndname(ielem, ibin))//'TP', bndiags(:, :, ibin, ielem, BNDIAGS_TP), pcols, lchnk)
              end if
            end if
          end if
        end if
      end do
    end do
    
    ! Output the particle diagnostics.
    do igroup = 1, NGROUP 
      call CARMAGROUP_Get(carma, igroup, rc, shortname=sname, is_cloud=is_cloud, do_drydep=do_drydep, ienconc=ienconc)
      if (rc < 0) call endrun('carma_output_diagnostics::CARMAGROUP_Get failed.')
      
      ! Gridbox average
      call outfld(trim(sname)//'ND', gpdiags(:, :, igroup, GPDIAGS_ND), pcols, lchnk) 
      call outfld(trim(sname)//'AD', gpdiags(:, :, igroup, GPDIAGS_AD), pcols, lchnk) 
      call outfld(trim(sname)//'MD', gpdiags(:, :, igroup, GPDIAGS_MD), pcols, lchnk) 
      call outfld(trim(sname)//'RE', gpdiags(:, :, igroup, GPDIAGS_RE), pcols, lchnk) 
      call outfld(trim(sname)//'RM', gpdiags(:, :, igroup, GPDIAGS_RM), pcols, lchnk) 
      call outfld(trim(sname)//'JN', gpdiags(:, :, igroup, GPDIAGS_JN), pcols, lchnk) 
      call outfld(trim(sname)//'MR', gpdiags(:, :, igroup, GPDIAGS_MR), pcols, lchnk) 
      call outfld(trim(sname)//'EX', gpdiags(:, :, igroup, GPDIAGS_EX), pcols, lchnk) 
      call outfld(trim(sname)//'OD', gpdiags(:, :, igroup, GPDIAGS_OD), pcols, lchnk) 
      call outfld(trim(sname)//'PA', gpdiags(:, :, igroup, GPDIAGS_PA), pcols, lchnk) 
      call outfld(trim(sname)//'AR', gpdiags(:, :, igroup, GPDIAGS_AR), pcols, lchnk) 
      call outfld(trim(sname)//'VM', gpdiags(:, :, igroup, GPDIAGS_VM), pcols, lchnk) 
      
      if (do_drydep) then
        do ibin = 1, NBIN
          call outfld(trim(btndname(igroup, ibin))//'VD', sbdiags(:, ibin, ienconc, SBDIAGS_VD), pcols, lchnk)
        end do
      end if
    end do
    
    ! Output the gas tendencies.
    do igas = 1, NGAS
      icnst = icnst4gas(igas)
  
      call outfld(gtndname(igas), ptend%q(:, :, icnst), pcols, lchnk) 
      
      ! Output the supersaturations.
      call outfld(trim(cnst_name(icnst))//'SI', gsdiags(:, :, igas, GSDIAGS_SI), pcols, lchnk) 
      call outfld(trim(cnst_name(icnst))//'SL', gsdiags(:, :, igas, GSDIAGS_SL), pcols, lchnk) 
      call outfld(trim(cnst_name(icnst))//'EI', gsdiags(:, :, igas, GSDIAGS_EI), pcols, lchnk) 
      call outfld(trim(cnst_name(icnst))//'EL', gsdiags(:, :, igas, GSDIAGS_EL), pcols, lchnk) 
      call outfld(trim(cnst_name(icnst))//'WT', gsdiags(:, :, igas, GSDIAGS_WT), pcols, lchnk) 
    end do
    
    ! Output the temperature tendency.
    if (carma_do_thermo) then
       call outfld('CRTT', ptend%s(:, :) / cpair, pcols, lchnk) 
    end if
              
    return
  end subroutine carma_output_diagnostics
       
    
  !! Calculate the emissions for CARMA aerosols. This is taken from
  !! the routine aerosol_emis_intr in aerosol_intr.F90 and dust_emis_intr in
  !! dust_intr.F90 by Phil Rasch.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine carma_emission_tend (state, ptend, cam_in, dt)
    use cam_history,   only: outfld
    use camsrfexch,       only: cam_in_t

    implicit none
    
    type(physics_state), intent(in )    :: state                !! physics state
    type(physics_ptend), intent(inout)  :: ptend                !! physics state tendencies
    type(cam_in_t),      intent(inout)  :: cam_in               !! surface inputs
    real(r8),            intent(in)     :: dt                   !! time step (s)

    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: igroup                  ! group index
    integer      :: ielem                   ! element index
    integer      :: ibin                    ! bin index
    integer      :: icnst                   ! consituent index
    real(r8)     :: tendency(pcols, pver)   ! constituent tendency (kg/kg/s)
    real(r8)     :: surfaceFlux(pcols)      ! constituent surface flux (kg/m^2/s)
    integer      :: cnsttype                ! constituent type
    integer      :: maxbin                  ! last prognostic bin
    integer      :: rc                      ! CARMA return code

    ! Initialize the return code.
    rc = 0

    ! Initialize the output tendency structure.
    call physics_ptend_init(ptend,state%psetcols, 'CARMA (emission)', lq=lq_carma)

    if (.not. carma_flag) return
    if (.not. carma_do_emission) return

    ncol = state%ncol
    lchnk = state%lchnk
    
    ! Provide emissions rates for particles.
    !
    ! NOTE: This can only be done for prognostic groups.
    do ielem = 1, NELEM
      call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
      if (rc < 0) call endrun('carma_drydep_tend::CARMAELEMENT_Get failed.')
      
      call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, maxbin=maxbin)
      if (rc < 0) call endrun('carma_drydep_tend::CARMAGROUP_Get failed.')
      
      if (cnsttype == I_CNSTTYPE_PROGNOSTIC) then
      
        do ibin = 1, NBIN

          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the paerent model.
          if (ibin <= maxbin) then
    
            icnst = icnst4elem(ielem, ibin)
          
            call CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
            if (rc < 0) call endrun('carma_emission_tend::CARMA_EmitParticle failed.')
          
            ! Add any surface flux here.
            cam_in%cflx(:ncol, icnst) = surfaceFlux(:ncol)
            call outfld(trim(cnst_name(icnst))//'SF', cam_in%cflx(:ncol, icnst), ncol, lchnk)
            
            ! For emissions into the atmosphere, put the emission here.
            ptend%q(:ncol, :pver, icnst) = tendency(:ncol, :pver)
            call outfld(trim(cnst_name(icnst))//'EM', ptend%q(:ncol, :, icnst), ncol, lchnk)
          end if
        enddo
      end if
    enddo
    
    ! No emissions rate is set up for gases, but it could be added here.

    return
  end subroutine carma_emission_tend 


  !! Calculate the wet deposition for the CARMA aerosols. This is taken from
  !! the routine aerosol_wet_int in aerosol_intr.F90 and dust_wet_intr in
  !! dust_intr.F90 by Phil Rasch.
  !! 
  !! Method: 
  !!  Use a modified version of the scavenging parameterization described in
  !!     Barth et al, 2000, JGR (sulfur cycle paper)
  !!     Rasch et al, 2001, JGR (INDOEX paper)
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine carma_wetdep_tend(state, ptend, dt,  pbuf, dlf, cam_out)
    use cam_history,   only: outfld
    use phys_control,  only: cam_physpkg_is
    use phys_grid,     only: get_lat_all_p, get_lon_all_p, get_rlat_all_p
    use wetdep,        only: clddiag, wetdepa_v1, wetdepa_v2
    use camsrfexch,       only: cam_out_t
    use physconst,     only: gravit
    
    implicit none

    real(r8),             intent(in)    :: dt              !! time step (s)
    type(physics_state),  intent(in )   :: state           !! physics state
    type(physics_ptend),  intent(inout) :: ptend           !! physics state tendencies
    type(physics_buffer_desc), pointer  :: pbuf(:)         !! physics buffer
    real(r8), intent(in)                :: dlf(pcols,pver) !! Detrainment of convective condensate
    type(cam_out_t),      intent(inout) :: cam_out         !! cam output to surface models

    ! local vars
    real(r8)                            :: rainmr(pcols,pver)     ! mixing ratio of rain within cloud volume
    real(r8)                            :: cldv(pcols,pver)       ! cloudy volume undergoing wet chem and scavenging
    real(r8)                            :: cldvcu(pcols,pver)   ! Convective precipitation area, top interface of current layer
    real(r8)                            :: cldvst(pcols,pver)   ! Stratiform precipitation area, top interface of current layer 
    integer                             :: ielem                  ! element index
    integer                             :: igroup                 ! group index
    integer                             :: ibin                   ! bin index
    integer                             :: icnst                  ! constituent index
    integer                             :: lat(pcols)             ! latitude indices
    real(r8)                            :: clat(pcols)            ! latitudes
    integer                             :: lon(pcols)             ! longtitude indices
    real(r8)                            :: conicw(pcols,pver)     ! convective in-cloud water
    real(r8)                            :: cmfdqr(pcols,pver)     ! convective production of rain
    real(r8)                            :: cldc(pcols,pver)       ! convective cloud fraction, currently empty
    real(r8)                            :: clds(pcols,pver)       ! Stratiform cloud fraction
    real(r8)                            :: evapc(pcols,pver)      ! Evaporation rate of convective precipitation
    real(r8)                            :: iscavt(pcols, pver)
    real(r8)                            :: scavt(pcols, pver)
    integer                             :: ixcldliq
    integer                             :: ixcldice
    real(r8)                            :: totcond(pcols, pver)   ! total condensate
    real(r8)                            :: solfac                 ! solubility factor
    real(r8)                            :: scavcoef               ! scavenging Coefficient
    logical                             :: do_wetdep
    integer                             :: ncol                   ! number of columns
    integer                             :: lchnk                  ! chunk identifier
    integer                             :: rc                     ! CARMA return code
    real(r8)                            :: z_scavcoef(pcols,pver) ! Dana and Hales coefficient (/mm)
    integer                             :: cnsttype               ! constituent type
    integer                             :: k
    real(r8)                            :: sflx(pcols)            ! Surface Flux (kg/m2/s)
    integer                             :: maxbin

    ! physics buffer 
    integer itim_old, ifld
    real(r8), pointer, dimension(:,:)   :: cldn                   ! cloud fraction
    real(r8), pointer, dimension(:,:)   :: cme
    real(r8), pointer, dimension(:,:)   :: prain
    real(r8), pointer, dimension(:,:)   :: evapr
    real(r8), pointer, dimension(:,:)   :: icwmrdp                ! in cloud water mixing ratio, deep convection
    real(r8), pointer, dimension(:,:)   :: rprddp                 ! rain production, deep convection
    real(r8), pointer, dimension(:,:)   :: icwmrsh                ! in cloud water mixing ratio, deep convection
    real(r8), pointer, dimension(:,:)   :: rprdsh                 ! rain production, deep convection
    real(r8), pointer, dimension(:,:,:) :: fracis                 ! fraction of transported species that are insoluble
    real(r8), pointer, dimension(:,:) ::  sh_frac  ! Shallow convective cloud fraction
    real(r8), pointer, dimension(:,:) ::  dp_frac  ! Deep convective cloud fraction
    real(r8), pointer, dimension(:,:) ::  evapcsh  ! Evaporation rate of shallow convective precipitation >=0.
    real(r8), pointer, dimension(:,:) ::  evapcdp  ! Evaporation rate of deep    convective precipitation >=0.

    ! Initialize the return code.
    rc = 0
    
    ! Initialize the output tendency structure.
    call physics_ptend_init(ptend,state%psetcols, 'CARMA (wetdep)', lq=lq_carma)

    if (.not. carma_flag) return
    if (.not. carma_do_wetdep) return    

    ncol = state%ncol
    lchnk = state%lchnk

    call get_lat_all_p(lchnk, ncol, lat)
    call get_lon_all_p(lchnk, ncol, lon)
    call get_rlat_all_p(lchnk, ncol, clat)

    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()
    
    call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cldn, (/1,1,itim_old/),(/pcols,pver,1/))
    call pbuf_get_field(pbuf, pbuf_get_index('QME'), cme )
    call pbuf_get_field(pbuf, pbuf_get_index('PRAIN'), prain )
    call pbuf_get_field(pbuf, pbuf_get_index('NEVAPR'), evapr )
    call pbuf_get_field(pbuf, pbuf_get_index('FRACIS'), fracis )
    call pbuf_get_field(pbuf, pbuf_get_index('ICWMRDP'), icwmrdp )
    call pbuf_get_field(pbuf, pbuf_get_index('RPRDDP'), rprddp )
    call pbuf_get_field(pbuf, pbuf_get_index('ICWMRSH'), icwmrsh )
    call pbuf_get_field(pbuf, pbuf_get_index('RPRDSH'), rprdsh )

    ! sum deep and shallow convection contributions
    conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
    cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)

    call pbuf_get_field(pbuf, pbuf_get_index('SH_FRAC'), sh_frac )
    call pbuf_get_field(pbuf, pbuf_get_index('DP_FRAC'), dp_frac )
    call pbuf_get_field(pbuf, pbuf_get_index('NEVAPR_SHCU'), evapcsh )
    call pbuf_get_field(pbuf, pbuf_get_index('NEVAPR_DPCU'), evapcdp )
    
    cldc(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:) ! Sungsu included this.
    evapc(:ncol,:) = evapcsh(:ncol,:) + evapcdp(:ncol,:) ! Sungsu included this.
    clds(:ncol,:)  = cldn(:ncol,:) - cldc(:ncol,:)       ! Stratiform cloud fraction


   cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)
   
    !   fields needed for wet scavenging
    call clddiag( state%t, state%pmid, state%pdel, cmfdqr, evapc, cldn, cldc, clds, cme, evapr, prain, &
         cldv, cldvcu, cldvst, rainmr, ncol )

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + &
         state%q(:ncol,:,ixcldice)

    ! Iterate over each particle and calculate a tendency from wet
    ! scavenging for it.
    do ielem = 1, NELEM
    
      ! NOTE: This can only be done for prognistic groups.
    
      call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
      if (rc < 0) call endrun('carma_wetdep_tend::CARMAELEMENT_Get failed.')
      
      call CARMAGROUP_Get(carma, igroup, rc, cnsttype=cnsttype, do_wetdep=do_wetdep, &
        solfac=solfac, scavcoef=scavcoef, maxbin=maxbin)
      if (rc < 0) call endrun('carma_wetdep_tend::CARMAGROUP_Get failed.')
      
      if ((do_wetdep) .and. (cnsttype == I_CNSTTYPE_PROGNOSTIC)) then
      
        do ibin = 1, NBIN
        
          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the parent model.
          if (ibin <= maxbin) then
    
            icnst = icnst4elem(ielem, ibin)
            
            scavt    = 0._r8
            
            ! The scavenging coefficient might be calculated as a function of
            ! the aerosol bin at each grid point. However, for now, we will just
            ! use a constant value for each group.
            z_scavcoef(:, :) = scavcoef
    
            if (cam_physpkg_is('cam5')) then

              call wetdepa_v2(state%t, &
                           state%pmid, &
                           state%q, &
                           state%pdel, &
                           cldn, &
                           cldc, &
                           cmfdqr, &
                           evapc, &
                           conicw, &
                           prain, &
                           cme, &
                           evapr, &
                           totcond, &
                           state%q(:, :, icnst), & 
                           dt, &
                           scavt, &
                           iscavt, &
                           cldv, &
                           cldvcu, &
                           cldvst, &
                           dlf, & 
                           fracis(:, :, icnst), & 
                           solfac, &
                           ncol, &
                           z_scavcoef)
                          
            else if (cam_physpkg_is('cam4')) then
 
              call wetdepa_v1(state%t, &
                           state%pmid, &
                           state%q, &
                           state%pdel, &
                           cldn, &
                           cldc, &
                           cmfdqr, &
                           conicw, &
                           prain, &
                           cme, &
                           evapr, &
                           totcond, &
                           state%q(:, :, icnst), & 
                           dt, &
                           scavt, &
                           iscavt, &
                           cldv, &
                           fracis(:, :, icnst), & 
                           solfac, &
                           ncol, &
                           z_scavcoef)
            else
            
              call endrun('carma_wetdep_tend:: No wet deposition routine is available for this configuration.')
            end if
                         
            ptend%q(:, :, icnst) = scavt
            call outfld(trim(cnst_name(icnst))//'WD', ptend%q(:, :, icnst), pcols, lchnk)

	          !
	          ! ptend%q(kg/kg air/s) * pdel(Pa) / gravit (m/s2) => (kg/m2/s)
	          !  note: 1Pa = 1 kg air * (m/s2) / m2
            sflx(:) = 0._r8
            
            do k = 1,pver
              sflx(:ncol) = sflx(:ncol) - ptend%q(:ncol, k, icnst) * state%pdel(:ncol,k) / gravit
            enddo

            call outfld(trim(cnst_name(icnst))//'SW', sflx, pcols, lchnk)

            ! Add this to the surface amount of the constituent
            call CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
 
          end if
        end do
      end if
    end do

    return
  end subroutine carma_wetdep_tend
  
  
  !! This routine creates files containing optical properties for each radiatively
  !! active particle type. These optical properties are used by the RRTMG radiation
  !! code to include the impact of CARMA particles in the radiative transfer
  !! calculation.
  !!
  !! NOTE: The format of this file is determined by the needs of the radiative tranfer
  !! code, so ideally a routine would exist in that module that could create a file
  !! with the proper format. Since that doesn't exist, we do it all here.
  subroutine CARMA_CreateOpticsFile(carma, rc)
    use radconstants, only : nswbands, nlwbands
    use wrap_nf
    use wetr, only         : getwetr
    
    implicit none

    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(out)                :: rc            !! return code, negative indicates failure

    ! Local variables
    integer                             :: igroup, ibin, iwave, irh
    integer                             :: irhswell
    integer                             :: ienconc
    real(kind=f)                        :: rho(NBIN), rhopwet
    real(kind=f)                        :: r(NBIN), rmass(NBIN), rlow(NBIN), rup(NBIN)
    real(kind=f)                        :: wave(NWAVE)
    complex(kind=f)                     :: refidx(NWAVE)
    character(len=CARMA_NAME_LEN)       :: name
    character(len=CARMA_SHORT_NAME_LEN) :: shortname
    logical                             :: do_mie
    integer                             :: fid
    integer                             :: rhdim, lwdim, swdim
    integer                             :: rhvar, lwvar, swvar
    integer                             :: abs_lw_var
    integer                             :: ext_sw_var, ssa_sw_var, asm_sw_var
    integer                             :: omdim, andim, namedim 
    integer                             :: omvar, anvar, namevar 
    integer                             :: dimids(2)
    integer                             :: denvar, slogvar, dryrvar, rminvar, rmaxvar, hygrovar, ntmvar
    real(kind=f)                        :: abs_lw(NMIE_RH, nlwbands)
    real(kind=f)                        :: ext_sw(NMIE_RH, nswbands)
    real(kind=f)                        :: ssa_sw(NMIE_RH, nswbands)
    real(kind=f)                        :: asm_sw(NMIE_RH, nswbands)
    character(len=8)                    :: c_name                   ! constituent name  
    character(len=32)                   :: aer_name                 ! long enough for both aername and name  
    character(len=255)                  :: filepath
    real(kind=f)                        :: rwet
    real(kind=f)                        :: Qext
    real(kind=f)                        :: Qsca
    real(kind=f)                        :: asym
    integer                             :: start_text(2), count_text(2)
    integer                             :: sw_r_refidx_var, sw_i_refidx_var, lw_r_refidx_var, lw_i_refidx_var
    integer                             :: nrh
    integer                             :: cnsttype               ! constituent type
    integer                             :: maxbin                 ! last prognostic bin
    integer                             :: LUNOPRT              ! logical unit number for output
    logical                             :: do_print             ! do print output?
    integer                             :: ret
    
    
    ! Assume success.
    rc = 0
    
    ! Get the wavelength structure.
    call CARMA_GET(carma, rc, wave=wave, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun('carma_CreateOpticsFile::CARMA_Get failed.')
    
    ! Process each group that is defined in the model.
    do igroup = 1, NGROUP
    
      ! Get the necessary group properties.
      call CARMAGROUP_Get(carma, igroup, rc, do_mie=do_mie, name=name, shortname=shortname, r=r, &
                          rlow=rlow, rup=rup, rmass=rmass, refidx=refidx, irhswell=irhswell, &
                          ienconc=ienconc, cnsttype=cnsttype, maxbin=maxbin)
      if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAGROUP_Get failed.')
      
      ! Are we supposed to do the mie calculation for this group?
      if ((do_mie) .and. (cnsttype == I_CNSTTYPE_PROGNOSTIC)) then
      
        call CARMAELEMENT_Get(carma, ienconc, rc, rho=rho)
        if (rc < 0) call endrun('carma_CreateOpticsFile::CARMAELEMENT_Get failed.')
      
        ! A file needs to be created for each bin.
        do ibin = 1, NBIN
        
          ! Bins past maxbin are treated as diagnostic even if the group
          ! is prognostic and thus are not advected in the paerent model.
          if (ibin <= maxbin) then
    
            write(c_name, '(A, I2.2)') trim(shortname), ibin
          
            ! Construct the path to the file. Each model will have its own subdirectory
            ! where the optical property files are stored.
            filepath = trim(carma_model) // '_' // trim(c_name) // '_rrtmg.nc'
            
            if (do_print) write(LUNOPRT,*) 'Creating CARMA optics file ... ', trim(filepath)
    
            ! Create the file.
            call wrap_create(filepath, NF90_CLOBBER, fid)
            
            ! For non-hygroscopic, only use 1 RH value.
            if (irhswell /= 0) then
              nrh = NMIE_RH
            else
              nrh = min(NMIE_RH, 1)
            end if
              
            ! Define the dimensions: rh, lwbands, swbands
            call wrap_def_dim(fid, 'rh_idx',  nrh,  rhdim)
            call wrap_def_dim(fid, 'lw_band', nlwbands, lwdim)
            call wrap_def_dim(fid, 'sw_band', nswbands, swdim)

            write(LUNOPRT,*) "Defined rh_idx, lw_band, and sw_band dims."

            dimids(1) = rhdim
            call wrap_def_var(fid, 'rh',  NF90_DOUBLE, 1, dimids(1:1), rhvar)
            
            dimids(1) = lwdim
            call wrap_def_var(fid, 'lw_band', NF90_DOUBLE, 1, dimids(1:1), lwvar)
            
            dimids(1) = swdim
            call wrap_def_var(fid, 'sw_band', NF90_DOUBLE, 1, dimids(1:1), swvar)

            write(LUNOPRT,*) "Defined rh_idx, lw_band, and sw_band vars."
            
            call wrap_put_att_text(fid, rhvar, 'units', 'fraction') 
            call wrap_put_att_text(fid, lwvar, 'units', 'm') 
            call wrap_put_att_text(fid, swvar, 'units', 'm') 
  
            call wrap_put_att_text(fid, rhvar, 'long_name', 'relative humidity')
            call wrap_put_att_text(fid, lwvar, 'long_name', 'longwave bands')
            call wrap_put_att_text(fid, swvar, 'long_name', 'shortwave bands')
            
            ! Define the variables: abs_lw, ext_sw, ssa_sw, asm_sw
            dimids(1) = rhdim
            dimids(2) = lwdim
            call wrap_def_var(fid, 'abs_lw', NF90_DOUBLE, 2, dimids, abs_lw_var)

            write(LUNOPRT,*) "Defined abs_lw."

            call wrap_put_att_text(fid, abs_lw_var, 'units', 'meter^2 kilogram^-1') 
            
            dimids(1) = rhdim
            dimids(2) = swdim
            call wrap_def_var(fid, 'ext_sw', NF90_DOUBLE, 2, dimids, ext_sw_var)
            call wrap_def_var(fid, 'ssa_sw', NF90_DOUBLE, 2, dimids, ssa_sw_var)
            call wrap_def_var(fid, 'asm_sw', NF90_DOUBLE, 2, dimids, asm_sw_var)
            
            write(LUNOPRT,*) "Defined ext_sw, ssa_sw, and asm_sw."

            call wrap_put_att_text(fid, ssa_sw_var, 'units', 'fraction') 
            call wrap_put_att_text(fid, ext_sw_var, 'units', 'meter^2 kilogram^-1') 
            call wrap_put_att_text(fid, asm_sw_var, 'units', '-') 
            
            ! Define the variables for the refractive indicies.
            dimids(1) = swdim
            call wrap_def_var(fid, 'refindex_real_aer_sw', NF90_DOUBLE, 1, dimids(1:1), sw_r_refidx_var)
            call wrap_def_var(fid, 'refindex_im_aer_sw',   NF90_DOUBLE, 1, dimids(1:1), sw_i_refidx_var)
            
            write(LUNOPRT,*) "Defined lw refindex."

            dimids(1) = lwdim
            call wrap_def_var(fid, 'refindex_real_aer_lw', NF90_DOUBLE, 1, dimids(1:1), lw_r_refidx_var)
            call wrap_def_var(fid, 'refindex_im_aer_lw',   NF90_DOUBLE, 1, dimids(1:1), lw_i_refidx_var)
  
            write(LUNOPRT,*) "Defined sw refindex."

            call wrap_put_att_text(fid, sw_r_refidx_var, 'units', '-') 
            call wrap_put_att_text(fid, sw_i_refidx_var, 'units', '-') 
            call wrap_put_att_text(fid, lw_r_refidx_var, 'units', '-') 
            call wrap_put_att_text(fid, lw_i_refidx_var, 'units', '-') 
  
            call wrap_put_att_text(fid, sw_r_refidx_var, 'long_name', 'real refractive index of aerosol - shortwave') 
            call wrap_put_att_text(fid, sw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - shortwave') 
            call wrap_put_att_text(fid, lw_r_refidx_var, 'long_name', 'real refractive index of aerosol - longwave') 
            call wrap_put_att_text(fid, lw_i_refidx_var, 'long_name', 'imaginary refractive index of aerosol - longwave') 
           
            
            ! Define fields that define the aerosol properties.
            call wrap_def_dim(fid, 'opticsmethod_len',  32, omdim)
            dimids(1) = omdim
            call wrap_def_var(fid, 'opticsmethod',  NF90_CHAR, 1, dimids(1:1), omvar)
  
            write(LUNOPRT,*) "Defined omdim."

            call wrap_def_dim(fid, 'namelength',  20, andim)
            dimids(1) = andim
            call wrap_def_var(fid, 'aername',  NF90_CHAR, 1, dimids(1:1), anvar)
  
            write(LUNOPRT,*) "Defined aername."

            call wrap_def_dim(fid, 'name_len',  32, namedim)
            dimids(1) = namedim
            call wrap_def_var(fid, 'name',  NF90_CHAR, 1, dimids(1:1), namevar)
  
            write(LUNOPRT,*) "Defined name."

            call wrap_def_var(fid, 'density',            NF90_DOUBLE, 0, dimids(1:0), denvar)
            call wrap_def_var(fid, 'sigma_logr',         NF90_DOUBLE, 0, dimids(1:0), slogvar)
            call wrap_def_var(fid, 'dryrad',             NF90_DOUBLE, 0, dimids(1:0), dryrvar)
            call wrap_def_var(fid, 'radmin_aer',         NF90_DOUBLE, 0, dimids(1:0), rminvar)
            call wrap_def_var(fid, 'radmax_aer',         NF90_DOUBLE, 0, dimids(1:0), rmaxvar)
            call wrap_def_var(fid, 'hygroscopicity',     NF90_DOUBLE, 0, dimids(1:0), hygrovar)
            call wrap_def_var(fid, 'num_to_mass_ratio',  NF90_DOUBLE, 0, dimids(1:0), ntmvar)
            
            call wrap_put_att_text(fid, denvar,   'units', 'kg m^-3') 
            call wrap_put_att_text(fid, slogvar,  'units', '-') 
            call wrap_put_att_text(fid, dryrvar,  'units', 'm') 
            call wrap_put_att_text(fid, rminvar,  'units', 'm') 
            call wrap_put_att_text(fid, rmaxvar,  'units', 'm') 
            call wrap_put_att_text(fid, hygrovar, 'units', '-') 
            call wrap_put_att_text(fid, ntmvar,   'units', 'kg^-1') 
  
            call wrap_put_att_text(fid, denvar,   'long_name', 'aerosol material density') 
            call wrap_put_att_text(fid, slogvar,  'long_name', 'geometric standard deviation of aerosol') 
            call wrap_put_att_text(fid, dryrvar,  'long_name', 'dry number mode radius of aerosol') 
            call wrap_put_att_text(fid, rminvar,  'long_name', 'minimum dry radius of aerosol for bin') 
            call wrap_put_att_text(fid, rmaxvar,  'long_name', 'maximum dry radius of aerosol for bin') 
            call wrap_put_att_text(fid, hygrovar, 'long_name', 'hygroscopicity of aerosol') 
            call wrap_put_att_text(fid, ntmvar,   'long_name', 'ratio of number to mass of aerosol') 
            
  
            write(LUNOPRT,*) "Defined all variables."

            ! End the defintion phase of the netcdf file.      
            call wrap_enddef(fid)
            
            
            ! Write out the dimensions.
            call wrap_put_var_realx(fid, rhvar, mie_rh(:nrh))
            call wrap_put_var_realx(fid, lwvar, wave(:nlwbands) * 1e-2_f)
            call wrap_put_var_realx(fid, swvar, wave(nlwbands+1:) * 1e-2_f)
            
            ! Write out the refractive indicies.
            call wrap_put_var_realx(fid, sw_r_refidx_var, real(refidx(nlwbands+1:)))
            call wrap_put_var_realx(fid, sw_i_refidx_var, aimag(refidx(nlwbands+1:)))
            call wrap_put_var_realx(fid, lw_r_refidx_var, real(refidx(:nlwbands)))
            call wrap_put_var_realx(fid, lw_i_refidx_var, aimag(refidx(:nlwbands)))
            
  
            ! Pad the names out with spaces.
            aer_name = '                                '
            aer_name(1:len(trim(c_name))) = c_name

            start_text(1) = 1
            count_text(1) = 32
            call wrap_put_vara_text(fid, namevar, start_text, count_text, (/ aer_name /))
            count_text(1) = 20
            call wrap_put_vara_text(fid, anvar, start_text, count_text, (/ aer_name /))
            
            ! These fields control whether the particle is treated as a CCN. For now,
            ! set these so that CARMA particles are not considered as CCN by the
            ! CAM microphysics.
            if (irhswell /= 0) then
              count_text(1) = len('hygroscopic                     ')
              call wrap_put_vara_text(fid, omvar, start_text, count_text, (/ 'hygroscopic                     ' /))
            else
              count_text(1) = len('insoluble                       ')
              call wrap_put_vara_text(fid, omvar, start_text, count_text, (/ 'insoluble                       ' /))
            end if
            
            call wrap_put_var_realx(fid, denvar,   (/ rho(ibin) * 1e-3_f / 1e-6_f /))
            call wrap_put_var_realx(fid, slogvar,  (/ 0._f /))
            call wrap_put_var_realx(fid, dryrvar,  (/ r(ibin) * 1e-2_f /))
            call wrap_put_var_realx(fid, rminvar,  (/ rlow(ibin) * 1e-2_f /))
            call wrap_put_var_realx(fid, rmaxvar,  (/ rup(ibin) * 1e-2_f /))
            call wrap_put_var_realx(fid, hygrovar, (/ 0._f /))
            call wrap_put_var_realx(fid, ntmvar,   (/ 1._f / rmass(ibin) / 1e-3_f /))
           
            ! Iterate over a range of relative humidities, since the particle may swell
            ! with relative humidity which will change its optical properties.
            do irh = 1, nrh
            
              ! Determine the wet radius.
              call getwetr(carma, igroup, mie_rh(irh), r(ibin), rwet, rho(ibin), rhopwet, rc)
              if (rc < 0) call endrun('carma_CreateOpticsFile::wetr failed.')
              
              ! Calculate at each wavelength.
              do iwave = 1, NWAVE
write(carma%f_LUNOPRT,*) "CARMA mie calc:  start  ", igroup, ibin, iwave, carma%f_wave(iwave), carma%f_group(igroup)%f_nmon(ibin)

        
                ! Using Mie code, calculate the optical properties: extinction coefficient,
                ! single scattering albedo and asymmetry factor.
                ! Assume the particle is homogeneous (no core).
                !
                ! NOTE: nmon, df, rmon and falpha are only used for fractal particles.
                call mie(carma, &
                         carma%f_group(igroup)%f_imiertn, &
                         rwet, &
                         carma%f_wave(iwave), &
                         carma%f_group(igroup)%f_nmon(ibin), &
                         carma%f_group(igroup)%f_df(ibin), &
                         carma%f_group(igroup)%f_rmon, &
                         carma%f_group(igroup)%f_falpha, &
                         carma%f_group(igroup)%f_refidx(iwave), &
                         Qext, &
                         Qsca, &
                         asym, &
                         rc)
                if (rc < 0) call endrun('carma_CreateOpticsFile::mie failed.')
write(carma%f_LUNOPRT,*) "CARMA mie calc:  done  ", Qext, Qsca, asym
    
              
                ! Calculate  the shortwave and longwave properties?
                !
                ! NOTE: miess is in cgs units, but the optics file needs to be in mks
                ! units, so perform the necessary conversions.
                if (iwave <= nlwbands) then
              
                  ! Longwave just needs absorption: abs_lw.
                  abs_lw(irh, iwave) = (Qext - Qsca) * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
                else
               
                  ! Shortwave needs extinction, single scattering albedo and asymmetry factor:
                  ! ext_sw, ssa_sw and asm_sw.
                  ext_sw(irh, iwave - nlwbands) = Qext * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
                  ssa_sw(irh, iwave - nlwbands) = Qsca / Qext
                  asm_sw(irh, iwave - nlwbands) = asym
                end if
              end do
            end do
            
            ! Write out the longwave fields.
            ret = nf90_put_var (fid, abs_lw_var, abs_lw(:nrh, :))
            if (ret/=NF90_NOERR) then
               write(iulog,*)'CARMA_CreateOpticsFile: error writing varid =', abs_lw_var
               call handle_error (ret)
            end if
            
            ! Write out the shortwave fields.
            ret = nf90_put_var (fid, ext_sw_var, ext_sw(:nrh, :))
            if (ret/=NF90_NOERR) then
               write(iulog,*)'CARMA_CreateOpticsFile: error writing varid =', ext_sw_var
               call handle_error (ret)
            end if
            ret = nf90_put_var (fid, ssa_sw_var, ssa_sw(:nrh, :))
            if (ret/=NF90_NOERR) then
               write(iulog,*)'CARMA_CreateOpticsFile: error writing varid =', ssa_sw_var
               call handle_error (ret)
            end if
            ret = nf90_put_var (fid, asm_sw_var, asm_sw(:nrh, :))
            if (ret/=NF90_NOERR) then
               write(iulog,*)'CARMA_CreateOpticsFile: error writing varid =', asm_sw_var
               call handle_error (ret)
            end if
            
            ! Close the file.
            call wrap_close(fid)
          end if  
        end do
      end if
    end do
    
    return
  end subroutine CARMA_CreateOpticsFile 
  
  
  !! This routine creates a file containing a reference temperature profile
  !! for use with fixed initialization.
  subroutine CARMA_CreateRefTFile(carma, filepath, lev, reft, rc, refh2o, refh2so4)
    use wrap_nf
    
    implicit none

    type(carma_type), intent(inout)     :: carma          !! the carma object
    character(len=*), intent(in)        :: filepath       !! the file path
    real(kind=f), intent(in)            :: lev(pver)      !! pressure levels
    real(kind=f), intent(in)            :: reft(pver)     !! reference temperature
    integer, intent(out)                :: rc             !! return code, negative indicates failure
    real(kind=f), optional, intent(in)  :: refh2o(pver)   !! reference water vapor
    real(kind=f), optional, intent(in)  :: refh2so4(pver) !! reference sulfuric acid

    ! Local variables
    integer                             :: fid
    integer                             :: levdim
    integer                             :: levvar, tvar, h2ovar, h2so4var
    integer                             :: dimids(2)
    
    
    ! Assume success.
    rc = 0
    
    ! Create the file.
    call wrap_create(filepath, NF90_CLOBBER, fid)
    
    
    ! Define the dimensions: lev
    call wrap_def_dim(fid, 'lev',  pver,  levdim)
            
    dimids(1) = levdim
    call wrap_def_var(fid, 'lev',  NF90_DOUBLE, 1, dimids(1:1), levvar)

    call wrap_put_att_text(fid, levvar, 'units', 'level') 
    call wrap_put_att_text(fid, levvar, 'long_name', 'hybrid level at midpoints (1000*(A+B))') 
    call wrap_put_att_text(fid, levvar, 'positive', 'down') 
    call wrap_put_att_text(fid, levvar, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate') 
    call wrap_put_att_text(fid, levvar, 'formula_terms', 'a: hyam b: hybm p0: P0 ps: PS') 
            
    ! Define the variables: T
    call wrap_def_var(fid, 'T', NF90_DOUBLE, 1, dimids(1:1), tvar)
            
    call wrap_put_att_text(fid, tvar, 'units', 'K') 
    call wrap_put_att_text(fid, tvar, 'long_name', 'Temperature')
    
    if ((carma%f_igash2o /= 0) .and. present(refh2o)) then
      call wrap_def_var(fid, 'Q', NF90_DOUBLE, 1, dimids(1:1), h2ovar)
            
      call wrap_put_att_text(fid, h2ovar, 'units', 'kg/kg') 
      call wrap_put_att_text(fid, h2ovar, 'long_name', 'Specific Humidity')
    end if

    if ((carma%f_igash2so4 /= 0) .and. present(refh2so4)) then
      call wrap_def_var(fid, 'H2SO4', NF90_DOUBLE, 1, dimids(1:1), h2so4var)
            
      call wrap_put_att_text(fid, h2so4var, 'units', 'kg/kg') 
      call wrap_put_att_text(fid, h2so4var, 'long_name', 'H2SO4')
    end if
  
    ! End the defintion phase of the netcdf file.      
    call wrap_enddef(fid)
    
    
    ! Write out the dimensions.
    call wrap_put_var_realx(fid, levvar, lev)
    
    ! Write out the variables.
    call wrap_put_var_realx(fid, tvar, reft)

    if ((carma%f_igash2o /= 0) .and. present(refh2o)) then
      call wrap_put_var_realx(fid, h2ovar, refh2o(:))
    end if

    if ((carma%f_igash2so4 /= 0) .and. present(refh2so4)) then
      call wrap_put_var_realx(fid, h2so4var, refh2so4(:))
    end if
    
    ! Close the file.
    call wrap_close(fid)
    
    return
  end subroutine CARMA_CreateRefTFile
  
  
  !! Calculate the aerodynamic resistance for dry deposition.
  !!
  !! This is based upon Seinfeld and Pandis (1998) page 963, and
  !! is similar to the calcram routine in drydep_mod.F90;
  !! however, it enables independent determination of the aerodynamic
  !! resistance each surface type (land, ocean and ice).
  !!
  !! @author  Tianyi Fan
  !! @version Aug 2011
  subroutine CARMA_calcram(ustar, z0, pdel, pmid, tmid, obklen, ram)
    use shr_const_mod, only: shr_const_karman  
    use physconst,     only: rair, gravit

    implicit none 

    ! input and output argument
    real(r8), intent(in)    :: ustar     ! friction velocity
    real(r8), intent(in)    :: z0        ! roughness length
    real(r8), intent(in)    :: pdel      ! layer depth in pressure  [Pa]
    real(r8), intent(in)    :: pmid      ! layer mid-point pressure [Pa]
    real(r8), intent(in)    :: tmid      ! layer mid-point temperature [K]
    real(r8), intent(in)    :: obklen    ! Monin-Obukhov length
    real(r8), intent(out)   :: ram       ! aerodynamic resistance
  
    ! local varibles
    real(r8)                :: z         ! half the layer height
    real(r8)                :: psi       ! stability parameter for z
    real(r8)                :: psi0      ! stability parameter for z0
    real(r8)                :: nu        ! temparory variable 
    real(r8)                :: nu0       ! temparory variable
    real(r8), parameter     :: xkar = shr_const_karman
  
    
    ! Use half the layer height like Ganzefeld and Lelieveld, 1995
    z = pdel * rair * tmid / pmid / gravit / 2._r8
    
    if (obklen .eq. 0._r8) then
      psi  = 0._r8
      psi0 = 0._r8
    else
      psi  = min(max(z  / obklen, -1._r8), 1._r8)
      psi0 = min(max(z0 / obklen, -1._r8), 1._r8) 
    endif
    
    ! Stable
    if (psi > 0._r8) then
      ram = 1._r8 / xkar / ustar * (log(z / z0) + 4.7_r8 * (psi - psi0))
    
    ! Unstable
    else if (psi < 0._r8) then
      nu  = (1._r8 - 15._r8 *psi)**(.25_r8)
      nu0 = (1._r8 - 15._r8 *psi0)**(.25_r8)

      if (ustar /= 0._r8) then
        ram = 1._r8 / xkar / ustar * (log(z / z0) + &
              log(((nu0**2 + 1._r8) * (nu0 + 1._r8)**2) / &
                  ((nu**2  + 1._r8) * (nu  + 1._r8)**2)) + &
                  2._r8 * (atan(nu) - atan(nu0)))
      else
        ram = 0._r8
      end if
    
    ! Neutral
    else
      ram = 1._r8 / xkar / ustar * log(z / z0)
    end if
    
    return 
  end subroutine CARMA_calcram 
end module carma_intr
