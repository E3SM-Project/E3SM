!! The CARMA module contains an interface to the Community Aerosol and Radiation
!! Model for Atmospheres (CARMA) bin microphysical model [Turco et al. 1979; 
!! Toon et al. 1988]. This implementation has been customized to work within
!! other model frameworks, so although it can be provided with an array of
!! columns, it does not do horizontal transport and just does independent 1-D
!! calculations upon each column.
!!
!! The typical usage for the CARMA and CARMASTATE objects within a model would be:
!!>
!!   ! This first section of code is done during the parent model's initialzation,
!!   ! and there should be a unique CARMA object created for each thread of
!!   ! execution.
!!
!!   ! Create the CARMA object.   
!!   call CARMA_Create(carma, ...)
!!
!!   ! Define the microphysical components. 
!!   call CARMAGROUP_Create(carma, ...)      ! One or more calls
!!
!!   call CARMAELEMENT_Create(carma, ...)  ! One or more calls
!!
!!   call CARMASOLUTE_Create(carma, ...)    ! Zero or more calls
!!
!!   call CARMAGAS_Create(carma, ...)          ! Zero or more calls
!!
!!   ! Define the relationships for the microphysical processes. 
!!   call CARMA_AddCoagulation(carma, ...)    ! Zero or more calls
!!   call CARMA_AddGrowth(carma, ...)         ! Zero or more calls
!!   call CARMA_AddNucleation(carma, ...)     ! Zero or more calls
!!
!!   ! Initialize things that are state and timestep independent.
!!   call CARMA_Initialize(carma, ...)
!!
!!   ...
!!   
!!   ! This section of code is within the parent model's timing loop.
!!   !
!!   ! NOTE: If using OPEN/MP, then each thread will execute one of
!!   ! of these loops per column of data. To avoid having to destroy
!!   ! the CARMASTATE object, a pool of CARMASTATE objects could be
!!   ! created so that there is one per thread and then the
!!   ! CARMA_Destroy() could be called after all columns have been
!!   ! processed.
!!
!!   ! Initialize CARMA for this model state and timestep.
!!   call CARMASTATE_Create(cstate, carma, ...)
!!
!!   ! Set the model state for each bin and gas.
!!   call CARMASTATE_SetBin(cstate, ...)          ! One call for each bin
!!   call CARMASTATE_SetGas(cstate, ...)          ! One call for each gas
!!
!!   ! Calculate the new state
!!   call CARMASTATE_Step(cstate, ...)
!!
!!   ! Get the results to return back to the parent model.
!!   call CARMASTATE_GetBin(cstate, ...)      ! One call for each Bin
!!   call CARMASTATE_GetGas(cstate, ...)      ! One call for each gas
!!   call CARMASTATE_GetState(cstate, ...)    ! Zero or one calls
!!
!!   ! (optional) Deallocate arrays that are not needed beyond this timestep.
!!   call CARMASTATE_Destroy(cstate)
!!
!!   ...
!!
!!   ! This section of code is done during the parent model's cleanup.
!!
!!   ! Deallocate all arrays.
!!   call CARMA_Destroy(carma)
!!<
!!
!!  @version Feb-2009 
!!  @author  Chuck Bardeen, Pete Colarco, Jamie Smith 
!
! NOTE: Documentation for this code can be generated automatically using f90doc,
! which is freely available from:
!   http://erikdemaine.org/software/f90doc/
! Comment lines with double comment characters are processed by f90doc, and there are
! some special characters added to the comments to control the documentation process.
! In addition to the special characters mentioned in the f990doc documentation, html
! formatting tags (e.g. <i></i>, <sup></sup>, ...) can also be added to the f90doc
! comments.
module carma_mod

  ! This module maps the parents models constants into the constants need by CARMA. NOTE: CARMA
  ! constants are in CGS units, while the parent models are typically in MKS units.
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMA_AddCoagulation
  public CARMA_AddGrowth
  public CARMA_AddNucleation
  public CARMA_Create
  public CARMA_Destroy
  public CARMA_Get
  public CARMA_Initialize

contains

  ! These are the methods that provide the interface between the parent model and the CARMA
  ! microphysical model. There are many other methods that are not in this file that are
  ! used to implement the microphysical calculations needed by the CARMA model. These other
  ! methods are in effect private methods of the CARMA module, but are in individual files
  ! since that is the way that CARMA has traditionally been structured and where users may
  ! want to extend or replace code to affect the microphysics.

  !! Creates the CARMA object and allocates arrays to store configuration information
  !! that will follow from the CARMA_AddXXX() methods. When the CARMA object is no longer
  !! needed, the CARMA_Destroy() method should be used to clean up any allocations
  !! that have happened. If LUNOPRT is specified, then the logical unit should be open and
  !! ready for output. The caller is responsible for closing the LUNOPRT logical unit
  !! after the CARMA object has been destroyed.
  !!
  !!  @version Feb-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
    LUNOPRT, wave, dwave, do_wave_emit)

    type(carma_type), intent(out)      :: carma     !! the carma object
    integer, intent(in)                :: NBIN      !! number of radius bins per group
    integer, intent(in)                :: NELEM     !! total number of elements
    integer, intent(in)                :: NGROUP    !! total number of groups
    integer, intent(in)                :: NSOLUTE   !! total number of solutes
    integer, intent(in)                :: NGAS      !! total number of gases
    integer, intent(in)                :: NWAVE     !! number of wavelengths
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    integer, intent(in), optional      :: LUNOPRT   !! logical unit number for output
    real(kind=f), intent(in), optional :: wave(NWAVE)  !! wavelength centers (cm)
    real(kind=f), intent(in), optional :: dwave(NWAVE) !! wavelength width (cm)
    logical, intent(in), optional      :: do_wave_emit(NWAVE) !! do emission in band?

    ! Local Varaibles      
    integer                            :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Save off the logic unit used for output if one was provided. If one was provided,
    ! then assume that CARMA can print output.
    if (present(LUNOPRT)) then 
      carma%f_LUNOPRT = LUNOPRT
      carma%f_do_print = .TRUE.
    end if
    
    ! Save the defintion of the number of comonents involved in the microphysics.
    carma%f_NGROUP  = NGROUP 
    carma%f_NELEM   = NELEM
    carma%f_NBIN    = NBIN
    carma%f_NGAS    = NGAS
    carma%f_NSOLUTE = NSOLUTE      
    carma%f_NWAVE   = NWAVE


    ! Allocate tables for the groups.
    allocate( &
      carma%f_group(NGROUP), &
      carma%f_icoag(NGROUP, NGROUP), &
      carma%f_inucgas(NGROUP), &
      stat=ier) 
    if(ier /= 0) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Create: ERROR allocating groups, NGROUP=", &
        carma%f_NGROUP, ", status=", ier
      rc = RC_ERROR
      return
    endif
    
    ! Initialize
    carma%f_icoag(:, :)  = 0
    carma%f_inucgas(:)   = 0    
    
    ! Allocate tables for the elements.
    allocate( &
      carma%f_element(NELEM), &
      carma%f_igrowgas(NELEM), &
      carma%f_inuc2elem(NELEM, NELEM), &
      carma%f_inucproc(NELEM, NELEM), &
      carma%f_ievp2elem(NELEM), &
      carma%f_nnuc2elem(NELEM), &
      carma%f_nnucelem(NELEM), &
      carma%f_inucelem(NELEM,NELEM*NGROUP), &
      carma%f_if_nuc(NELEM,NELEM), &
      carma%f_rlh_nuc(NELEM, NELEM), &
      carma%f_icoagelem(NELEM, NGROUP), &
      carma%f_icoagelem_cm(NELEM, NGROUP), &
      stat=ier) 
    if(ier /= 0) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Create: ERROR allocating elements, NELEM=", &
        carma%f_NELEM, ", status=", ier
      rc = RC_ERROR
      return
    endif
    
    ! Initialize
    carma%f_igrowgas(:) = 0
    carma%f_inuc2elem(:,:) = 0
    carma%f_inucproc(:,:) = 0
    carma%f_ievp2elem(:) = 0
    carma%f_nnuc2elem(:) = 0
    carma%f_nnucelem(:) = 0
    carma%f_inucelem(:,:) = 0
    carma%f_if_nuc(:,:) = .FALSE.
    carma%f_rlh_nuc(:,:) = 0._f
    carma%f_icoagelem(:,:) = 0
    carma%f_icoagelem_cm(:,:) = 0
    
    
    ! Allocate tables for the bins.
    allocate( &
      carma%f_inuc2bin(NBIN,NGROUP,NGROUP), &
      carma%f_ievp2bin(NBIN,NGROUP,NGROUP), &
      carma%f_nnucbin(NGROUP,NBIN,NGROUP), &
      carma%f_inucbin(NBIN*NGROUP,NGROUP,NBIN,NGROUP), &
      carma%f_diffmass(NBIN, NGROUP, NBIN, NGROUP), &
      carma%f_volx(NGROUP,NGROUP,NGROUP,NBIN,NBIN), &
      carma%f_ilow(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_jlow(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_iup(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_jup(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_npairl(NGROUP,NBIN), &
      carma%f_npairu(NGROUP,NBIN), &
      carma%f_iglow(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_jglow(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_igup(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_jgup(NGROUP,NBIN,NBIN*NBIN), &
      carma%f_kbin(NGROUP,NGROUP,NGROUP,NBIN,NBIN), &
      carma%f_pkernel(NBIN,NBIN,NGROUP,NGROUP,NGROUP,6), &
      carma%f_pratt(3,NBIN,NGROUP), &
      carma%f_prat(4,NBIN,NGROUP), &
      carma%f_pden1(NBIN,NGROUP), &
      carma%f_palr(4,NGROUP), &
      stat=ier) 
    if(ier /= 0) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Create: ERROR allocating bins, NBIN=", &
        carma%f_NBIN, ", status=", ier
      rc = RC_ERROR
      return
    endif
    
    ! Initialize
    carma%f_inuc2bin(:,:,:) = 0
    carma%f_ievp2bin(:,:,:) = 0
    carma%f_nnucbin(:,:,:) = 0
    carma%f_inucbin(:,:,:,:) = 0
    carma%f_diffmass(:, :, :, :) = 0._f
    carma%f_volx(:,:,:,:,:) = 0._f
    carma%f_ilow(:,:,:) = 0
    carma%f_jlow(:,:,:) = 0
    carma%f_iup(:,:,:) = 0
    carma%f_jup(:,:,:) = 0
    carma%f_npairl(:,:) = 0
    carma%f_npairu(:,:) = 0
    carma%f_iglow(:,:,:) = 0
    carma%f_jglow(:,:,:) = 0
    carma%f_igup(:,:,:) = 0
    carma%f_jgup(:,:,:) = 0
    carma%f_kbin(:,:,:,:,:) = 0._f
    carma%f_pkernel(:,:,:,:,:,:) = 0._f
    carma%f_pratt(:,:,:) = 0._f
    carma%f_prat(:,:,:) = 0._f
    carma%f_pden1(:,:) = 0._f
    carma%f_palr(:,:) = 0._f
      

    ! Allocate tables for solutes, if any are needed.
    if (NSOLUTE > 0) then
      allocate( &
        carma%f_solute(NSOLUTE), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Create: ERROR allocating solutes, NSOLUTE=", &
          carma%f_NSOLUTE, ", status=", ier 
        rc = RC_ERROR
        return
      endif
    end if
    
   
    ! Allocate tables for gases, if any are needed.
    if (NGAS > 0) then
      allocate( &
        carma%f_gas(NGAS), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Create: ERROR allocating gases, NGAS=", &
          carma%f_NGAS, ", status=", ier
        rc = RC_ERROR
        return
      endif
    end if
    
    
    ! Allocate tables for optical properties, if any are needed.
    if (NWAVE > 0) then
      allocate( &
        carma%f_wave(NWAVE), &
        carma%f_dwave(NWAVE), &
        carma%f_do_wave_emit(NWAVE), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Create: ERROR allocating wavelengths, NWAVE=", &
          carma%f_NWAVE, ", status=", ier
        rc = RC_ERROR
        return
      endif

      ! Initialize
      carma%f_do_wave_emit(:) = .TRUE.
      
      if (present(wave))  carma%f_wave(:)  = wave(:)
      if (present(dwave)) carma%f_dwave(:) = dwave(:)
      if (present(do_wave_emit)) carma%f_do_wave_emit(:) = do_wave_emit(:)
    end if
    
    return
 end subroutine CARMA_Create

  !! Called after the CARMA object has been created and the microphysics description has been
  !! configured. The optional flags control which microphysical processes are enabled and all of
  !! them default to FALSE. For a microphysical process to be active it must have been both
  !! configured (using a CARMA_AddXXX() method) and enabled here.
  !!
  !! NOTE: After initialization, the structure of the particle size bins is determined, and
  !! the resulting r, dr, rmass and dm can be retrieved with the CARMA_GetGroup() method.
  !!
  !!  @version Feb-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_Initialize(carma, rc, do_cnst_rlh, do_coag, do_detrain, do_fixedinit, &
      do_grow, do_incloud, do_explised, do_print_init, do_substep, do_thermo, do_vdiff, &
      do_vtran, do_drydep, vf_const, minsubsteps, maxsubsteps, maxretries, conmax, &
      do_pheat, do_pheatatm, dt_threshold, cstick, gsticki, gstickl, tstick, do_clearsky)
    type(carma_type), intent(inout)     :: carma         !! the carma object
    integer, intent(out)                :: rc            !! return code, negative indicates failure
    logical, intent(in), optional       :: do_cnst_rlh   !! use constant values for latent heats
                                                         !! (instead of varying with temperature)?
    logical, intent(in), optional       :: do_coag       !! do coagulation?
    logical, intent(in), optional       :: do_detrain    !! do detrainement?
    logical, intent(in), optional       :: do_fixedinit  !! do initialization from reference atm?
    logical, intent(in), optional       :: do_grow       !! do nucleation, growth and evaporation?
    logical, intent(in), optional       :: do_incloud    !! do incloud growth and coagulation?
    logical, intent(in), optional       :: do_explised   !! do sedimentation with substepping
    logical, intent(in), optional       :: do_substep    !! do substepping
    logical, intent(in), optional       :: do_print_init !! do prinit initializtion information
    logical, intent(in), optional       :: do_thermo     !! do thermodynamics
    logical, intent(in), optional       :: do_vdiff      !! do Brownian diffusion
    logical, intent(in), optional       :: do_vtran      !! do sedimentation
    logical, intent(in), optional       :: do_drydep     !! do dry deposition
    real(kind=f), intent(in), optional  :: vf_const      !! if specified and non-zero
                                                         !! constant fall velocity for all particles [cm/s]
    integer, intent(in), optional       :: minsubsteps   !! minimum number of substeps, default = 1
    integer, intent(in), optional       :: maxsubsteps   !! maximum number of substeps, default = 1
    integer, intent(in), optional       :: maxretries    !! maximum number of substep retries, default = 5
    real(kind=f), intent(in), optional  :: conmax        !! minimum relative concentration to consider, default = 1e-1
    logical, intent(in), optional       :: do_pheat      !! do particle heating
    logical, intent(in), optional       :: do_pheatatm   !! do particle heating of atmosphere
    real(kind=f), intent(in), optional  :: dt_threshold  !! convergence criteria for temperature [fraction]
    real(kind=f), intent(in), optional  :: cstick        !! accommodation coefficient - coagulation, default = 1.0
    real(kind=f), intent(in), optional  :: gsticki       !! accommodation coefficient - growth (ice), default = 0.93
    real(kind=f), intent(in), optional  :: gstickl       !! accommodation coefficient - growth (liquid), default = 1.0
    real(kind=f), intent(in), optional  :: tstick        !! accommodation coefficient - temperature, default = 1.0
    logical, intent(in), optional       :: do_clearsky   !! do clear sky growth and coagulation?
    
    ! Assume success.
    rc = RC_OK

    ! Set default values for control flags.
    carma%f_do_cnst_rlh   = .FALSE.
    carma%f_do_coag       = .FALSE.
    carma%f_do_detrain    = .FALSE.
    carma%f_do_fixedinit  = .FALSE.
    carma%f_do_grow       = .FALSE.
    carma%f_do_incloud    = .FALSE.
    carma%f_do_explised   = .FALSE.
    carma%f_do_pheat      = .FALSE.
    carma%f_do_pheatatm   = .FALSE.
    carma%f_do_print_init = .FALSE.
    carma%f_do_substep    = .FALSE.
    carma%f_do_thermo     = .FALSE.
    carma%f_do_vdiff      = .FALSE.
    carma%f_do_vtran      = .FALSE.
    carma%f_do_drydep     = .FALSE.
    carma%f_dt_threshold  = 0._f
    carma%f_cstick        = 1._f
    carma%f_gsticki       = 0.93_f
    carma%f_gstickl       = 1._f
    carma%f_tstick        = 1._f
    carma%f_do_clearsky   = .FALSE.

    ! Store off any control flag values that have been supplied.
    if (present(do_cnst_rlh))   carma%f_do_cnst_rlh   = do_cnst_rlh
    if (present(do_coag))       carma%f_do_coag       = do_coag
    if (present(do_detrain))    carma%f_do_detrain    = do_detrain
    if (present(do_fixedinit))  carma%f_do_fixedinit  = do_fixedinit
    if (present(do_grow))       carma%f_do_grow       = do_grow
    if (present(do_incloud))    carma%f_do_incloud    = do_incloud
    if (present(do_explised))   carma%f_do_explised   = do_explised
    if (present(do_pheat))      carma%f_do_pheat      = do_pheat
    if (present(do_pheatatm))   carma%f_do_pheatatm   = do_pheatatm
    if (present(do_print_init)) carma%f_do_print_init = (do_print_init .and. carma%f_do_print)
    if (present(do_substep))    carma%f_do_substep    = do_substep
    if (present(do_thermo))     carma%f_do_thermo     = do_thermo
    if (present(do_vdiff))      carma%f_do_vdiff      = do_vdiff
    if (present(do_vtran))      carma%f_do_vtran      = do_vtran 
    if (present(do_drydep))     carma%f_do_drydep     = do_drydep
    if (present(dt_threshold))  carma%f_dt_threshold  = dt_threshold
    if (present(cstick))        carma%f_cstick        = cstick
    if (present(gsticki))       carma%f_gsticki       = gsticki
    if (present(gstickl))       carma%f_gstickl       = gstickl
    if (present(tstick))        carma%f_tstick        = tstick
    if (present(do_clearsky))   carma%f_do_clearsky   = do_clearsky
    
    
    ! Setup the bin structure.
    call setupbins(carma, rc)
    if (rc < 0) return
    
    ! Substepping
    carma%f_minsubsteps = 1         ! minimum number of substeps
    carma%f_maxsubsteps = 1         ! maximum number of substeps
    carma%f_maxretries  = 1         ! maximum number of retries
    carma%f_conmax      = 1.e-1_f
    
    if (present(minsubsteps)) carma%f_minsubsteps = minsubsteps
    if (present(maxsubsteps)) carma%f_maxsubsteps = maxsubsteps
    if (present(maxretries))  carma%f_maxretries  = maxretries
    if (present(conmax))      carma%f_conmax      = conmax

    carma%f_do_step = .TRUE.
    
    ! Calculate the Optical Properties
    call CARMA_InitializeOptics(carma, rc)
    if (rc < 0) return

    ! If any of the processes have initialization that can be done without the state
    ! information, then perform that now. This will mostly be checking the configuration
    ! and setting up any tables based upon the configuration.
    if (carma%f_do_vtran .or. carma%f_do_coag)  then
      call CARMA_InitializeVertical(carma, rc, vf_const)
      if (rc < 0) return
    end if
    
    if (carma%f_do_coag) then
      call setupcoag(carma, rc)
      if (rc < 0) return
    end if
      
    if (carma%f_do_grow) then
      call CARMA_InitializeGrowth(carma, rc)
      if (rc < 0) return
    end if
      
    if (carma%f_do_thermo) then
      call CARMA_InitializeThermo(carma, rc)
       if (rc < 0) return
    end if
   
    return
  end subroutine CARMA_Initialize


  subroutine CARMA_InitializeGrowth(carma, rc)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)             :: rc
      
    ! Local Variables
    integer                            :: i
    logical                            :: bad_grid
    integer                            :: igroup   ! group index
    integer                            :: igas     ! gas index
    integer                            :: isol     ! solute index
    integer                            :: ielem    ! element index
    integer                            :: ibin     ! bin index
    integer                            :: igfrom
    integer                            :: igto
    integer                            :: ibto
    integer                            :: ieto
    integer                            :: ifrom
    integer                            :: iefrom
    integer                            :: jefrom
    integer                            :: ip
    integer                            :: jcore
    integer                            :: iecore
    integer                            :: im
    integer                            :: jnucelem
    integer                            :: inuc2
    integer                            :: neto
    integer                            :: jfrom
    integer                            :: j
    integer                            :: nnucb

    ! Define formats
    1 format(a,':  ',12i6)
    2 format(/,a,':  ',i6)
    3 format(a,a)
    4 format(a,':  ',1pe12.3)
    5 format(/,'Particle nucleation mapping arrays (setupnuc):')
    7 format(/,'Warning: nucleation cannot occur from group',i3, &
               '   bin',i3,'   into group',i3,'   (<inuc2bin> is zero)')
    
    
    ! Assume success.
    rc = RC_OK

    ! Compute radius-dependent terms used in PPM advection scheme
    do igroup = 1, carma%f_NGROUP
      do i = 2,carma%f_NBIN-1
        carma%f_pratt(1,i,igroup) = carma%f_group(igroup)%f_dm(i) / &
              ( carma%f_group(igroup)%f_dm(i-1) + carma%f_group(igroup)%f_dm(i) + carma%f_group(igroup)%f_dm(i+1) )
        carma%f_pratt(2,i,igroup) = ( 2._f*carma%f_group(igroup)%f_dm(i-1) + carma%f_group(igroup)%f_dm(i) ) / &
              ( carma%f_group(igroup)%f_dm(i+1) + carma%f_group(igroup)%f_dm(i) )
        carma%f_pratt(3,i,igroup) = ( 2._f*carma%f_group(igroup)%f_dm(i+1) + carma%f_group(igroup)%f_dm(i) ) / &
              ( carma%f_group(igroup)%f_dm(i-1) + carma%f_group(igroup)%f_dm(i) )
      enddo

      do i = 2,carma%f_NBIN-2
        carma%f_prat(1,i,igroup) = carma%f_group(igroup)%f_dm(i) / &
                ( carma%f_group(igroup)%f_dm(i) + carma%f_group(igroup)%f_dm(i+1) )
        carma%f_prat(2,i,igroup) = 2._f * carma%f_group(igroup)%f_dm(i+1) * carma%f_group(igroup)%f_dm(i) / &
               ( carma%f_group(igroup)%f_dm(i) + carma%f_group(igroup)%f_dm(i+1) )
        carma%f_prat(3,i,igroup) = ( carma%f_group(igroup)%f_dm(i-1) + carma%f_group(igroup)%f_dm(i) ) / &
               ( 2._f*carma%f_group(igroup)%f_dm(i) + carma%f_group(igroup)%f_dm(i+1) )
        carma%f_prat(4,i,igroup) = ( carma%f_group(igroup)%f_dm(i+2) + carma%f_group(igroup)%f_dm(i+1) ) / &
               ( 2._f*carma%f_group(igroup)%f_dm(i+1) + carma%f_group(igroup)%f_dm(i) )
        carma%f_pden1(i,igroup) = carma%f_group(igroup)%f_dm(i-1) + carma%f_group(igroup)%f_dm(i) + &
               carma%f_group(igroup)%f_dm(i+1) + carma%f_group(igroup)%f_dm(i+2)
      enddo

      if( carma%f_NBIN .gt. 1 )then
        carma%f_palr(1,igroup) = &
             (carma%f_group(igroup)%f_rmassup(1)-carma%f_group(igroup)%f_rmass(1)) / &
             (carma%f_group(igroup)%f_rmass(2)-carma%f_group(igroup)%f_rmass(1))
        carma%f_palr(2,igroup) = &
             (carma%f_group(igroup)%f_rmassup(1)/carma%f_group(igroup)%f_rmrat-carma%f_group(igroup)%f_rmass(1)) / &
             (carma%f_group(igroup)%f_rmass(2)-carma%f_group(igroup)%f_rmass(1))
        carma%f_palr(3,igroup) = &
             (carma%f_group(igroup)%f_rmassup(carma%f_NBIN-1)-carma%f_group(igroup)%f_rmass(carma%f_NBIN-1)) &
             / (carma%f_group(igroup)%f_rmass(carma%f_NBIN)-carma%f_group(igroup)%f_rmass(carma%f_NBIN-1))
        carma%f_palr(4,igroup) = &
             (carma%f_group(igroup)%f_rmassup(carma%f_NBIN)-carma%f_group(igroup)%f_rmass(carma%f_NBIN-1)) &
             / (carma%f_group(igroup)%f_rmass(carma%f_NBIN)-carma%f_group(igroup)%f_rmass(carma%f_NBIN-1))
      endif
    end do
    
    
    ! Check the nucleation mapping.
    !
    ! NOTE: This code was moved from setupnuc, because it is not dependent on the model's
    ! state. A small part of setupnuc which deals with scrit is state specific, and that was
    ! left in setupnuc.

    ! Bin mapping for nucleation : nucleation would transfer mass from particles
    ! in <ifrom,igfrom> into target bin <inuc2bin(ifrom,igfrom,igto)> in group
    ! <igto>.  The target bin is the smallest bin in the target size grid with
    ! mass exceeding that of nucleated particle.
    do igfrom = 1,carma%f_NGROUP    ! nucleation source group
      do igto = 1,carma%f_NGROUP        ! nucleation target group
        do ifrom = 1,carma%f_NBIN   ! nucleation source bin
  
          carma%f_inuc2bin(ifrom,igfrom,igto) = 0
  
          do ibto = carma%f_NBIN,1,-1        ! nucleation target bin
  
            if( carma%f_group(igto)%f_rmass(ibto) .ge. carma%f_group(igfrom)%f_rmass(ifrom) )then
              carma%f_inuc2bin(ifrom,igfrom,igto) = ibto
            endif
          enddo
        enddo
      enddo
    enddo

    ! Mappings for nucleation sources: 
    !
    !  <nnucelem(ielem)> is the number of particle elements that nucleate to
    !   particle element <ielem>.
    !
    !  <inuc2elem(jefrom,ielem)> are the particle elements that
    !   nucleate to particle element <ielem>, where 
    !   jefrom = 1,nnucelem(ielem).
    !
    !  <if_nuc(iefrom,ieto)> is true if nucleation transfers mass from element
    !   <iefrom> to element <ieto>.
    !
    !  <nnucbin(igfrom,ibin,igroup)> is the number of particle bins that nucleate
    !   to particles in bin <ibin,igroup> from group <igfrom>.
    !
    !  <inucbin(jfrom,igfrom,ibin,igto)> are the particle bins 
    !   that nucleate to particles in bin <ibin,igto>, where
    !   jfrom = 1,nnucbin(igfrom,ibin,igto).
    !
    !
    ! First, calculate <nnucelem(ielem)> and <if_nuc(iefrom,ieto)>
    ! based on <inucelem(jefrom,ielem)>
    do iefrom = 1,carma%f_NELEM
      do ieto = 1,carma%f_NELEM
        carma%f_if_nuc(iefrom,ieto) = .false.
      enddo
    enddo
    
    do ielem = 1,carma%f_NELEM
      carma%f_nnuc2elem(ielem) = 0
      
      do jefrom = 1,carma%f_NGROUP
        if( carma%f_inuc2elem(jefrom,ielem) .ne. 0 ) then
          carma%f_nnuc2elem(ielem) = carma%f_nnuc2elem(ielem) + 1
          carma%f_if_nuc(ielem,carma%f_inuc2elem(jefrom,ielem)) = .true.

      
          ! Also check for cases where neither the source or destinaton don't have cores (e.g.
          ! melting ice to water drops).
          if ((carma%f_group(carma%f_element(ielem)%f_igroup)%f_ncore .eq. 0) .and. &
              (carma%f_group(carma%f_element(carma%f_inuc2elem(jefrom,ielem))%f_igroup)%f_ncore .eq. 0)) then
      
            ! For particle concentration target elements, only count source elements
            ! that are also particle concentrations.
            carma%f_nnucelem(carma%f_inuc2elem(jefrom,ielem)) = carma%f_nnucelem(carma%f_inuc2elem(jefrom,ielem)) + 1
            carma%f_inucelem(carma%f_nnucelem(carma%f_inuc2elem(jefrom,ielem)),carma%f_inuc2elem(jefrom,ielem)) = ielem
          end if
        endif
      enddo
    enddo
    
    ! Next, enumerate and count elements that nucleate to cores.
    do igroup = 1,carma%f_NGROUP

      ip = carma%f_group(igroup)%f_ienconc    ! target particle number concentration element

      do jcore = 1,carma%f_group(igroup)%f_ncore

        iecore = carma%f_group(igroup)%f_icorelem(jcore)    ! target core element 
!        carma%f_nnucelem(iecore) = 0

        do iefrom = 1,carma%f_NELEM

          if( carma%f_if_nuc(iefrom,iecore) ) then
            carma%f_nnucelem(iecore) = carma%f_nnucelem(iecore) + 1
            carma%f_inucelem(carma%f_nnucelem(iecore),iecore) = iefrom
          endif
        enddo      ! iefrom=1,NELEM
      enddo        ! jcore=1,ncore
    enddo          ! igroup=1,NGROUP
    

    ! Now enumerate and count elements nucleating to particle concentration
    ! (itype=I_INVOLATILE and itype=I_VOLATILE) and core second moment
    ! (itype=I_COREMASS).  Elements with itype = I_VOLATILE are special because all
    ! nucleation sources for core elements in same group are also sources
    ! for the itype = I_VOLATILE element.
    do igroup = 1,carma%f_NGROUP
    
      ip = carma%f_group(igroup)%f_ienconc    ! target particle number concentration element
      im = carma%f_group(igroup)%f_imomelem   ! target core second moment element

!      carma%f_nnucelem(ip) = 0
!      if( im .ne. 0 )then
!        carma%f_nnucelem(im) = 0
!      endif

      do jcore = 1,carma%f_group(igroup)%f_ncore

        iecore = carma%f_group(igroup)%f_icorelem(jcore)       ! target core mass element

        do jnucelem = 1,carma%f_nnucelem(iecore)  ! elements nucleating to cores

          iefrom = carma%f_inucelem(jnucelem,iecore)  ! source
          
          ! For particle concentration target elements, only count source elements
          ! that are also particle concentrations.
          carma%f_nnucelem(ip) = carma%f_nnucelem(ip) + 1
          carma%f_inucelem(carma%f_nnucelem(ip),ip) = carma%f_group(carma%f_element(iefrom)%f_igroup)%f_ienconc

          if( im .ne. 0 )then
            carma%f_nnucelem(im) = carma%f_nnucelem(im) + 1
            carma%f_inucelem(carma%f_nnucelem(im),im) = iefrom
          endif
        enddo
      enddo       ! jcore=1,ncore
    enddo         ! igroup=1,NGROUP


    ! Now enumerate and count nucleating bins.
    do igroup = 1,carma%f_NGROUP    ! target group
      do ibin = 1,carma%f_NBIN    ! target bin
        do igfrom = 1,carma%f_NGROUP    ! source group

          carma%f_nnucbin(igfrom,ibin,igroup) = 0

          do ifrom = 1,carma%f_NBIN   ! source bin

            if( carma%f_inuc2bin(ifrom,igfrom,igroup) .eq. ibin ) then
              carma%f_nnucbin(igfrom,ibin,igroup) = carma%f_nnucbin(igfrom,ibin,igroup) + 1
              carma%f_inucbin(carma%f_nnucbin(igfrom,ibin,igroup),igfrom,ibin,igroup) = ifrom
            endif
          enddo
        enddo   ! igfrom=1,NGROUP
      enddo   ! ibin=1,NBIN=1,NGROUP
    enddo   ! igroup=1,NGROUP

    if (carma%f_do_print_init) then
      
      !  Report nucleation mapping arrays (should be 'write' stmts, of course)

      write(carma%f_LUNOPRT,*) ' '
      write(carma%f_LUNOPRT,*) 'Nucleation mapping arrays (setupnuc):'
      write(carma%f_LUNOPRT,*) ' '
      write(carma%f_LUNOPRT,*) 'Elements mapping:'
      
      do ielem = 1,carma%f_NELEM
        write(carma%f_LUNOPRT,*) 'ielem,nnucelem=',ielem,carma%f_nnucelem(ielem)
       
        if(carma%f_nnucelem(ielem) .gt. 0) then
          do jfrom = 1,carma%f_nnucelem(ielem)
            write(carma%f_LUNOPRT,*) 'jfrom,inucelem=  ',jfrom,carma%f_inucelem(jfrom,ielem)
          enddo
        endif
      enddo
      
      write(carma%f_LUNOPRT,*) ' '
      write(carma%f_LUNOPRT,*) 'Bin mapping:'
      
      do igfrom = 1,carma%f_NGROUP
        do igroup = 1,carma%f_NGROUP
          write(carma%f_LUNOPRT,*) ' '
          write(carma%f_LUNOPRT,*) 'Groups (from, to) = ', igfrom, igroup
          
          do ibin = 1,carma%f_NBIN
            nnucb = carma%f_nnucbin(igfrom,ibin,igroup)
            if(nnucb .eq. 0) write(carma%f_LUNOPRT,*) '  None for bin ',ibin
            if(nnucb .gt. 0) then
              write(carma%f_LUNOPRT,*) '  ibin,nnucbin=',ibin,nnucb
              write(carma%f_LUNOPRT,*) '   inucbin=',(carma%f_inucbin(j,igfrom,ibin,igroup),j=1,nnucb)
            endif
          enddo
        enddo
      enddo
    endif


    ! Check that values are valid.
    do ielem = 1, carma%f_NELEM

      if( carma%f_element(ielem)%f_isolute .gt. carma%f_NSOLUTE )then
        if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - component of isolute > NSOLUTE'
        rc = RC_ERROR
        return
      endif

      if( carma%f_ievp2elem(ielem) .gt. carma%f_NELEM )then
        if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - component of ievp2elem > NELEM'
        rc = RC_ERROR
        return
      endif

      ! Check that <isolute> is consistent with <ievp2elem>.
      if( carma%f_ievp2elem(ielem) .ne. 0 .and. carma%f_element(ielem)%f_itype .eq. I_COREMASS )then
        if( carma%f_element(ielem)%f_isolute .ne. carma%f_element(carma%f_ievp2elem(ielem))%f_isolute)then
          if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - isolute and ievp2elem are inconsistent'
          rc = RC_ERROR
          return
        endif
      endif

      ! Check that <isolute> is consistent with <inucgas>.
!      igas = carma%f_inucgas( carma%f_element(ielem)%f_igroup )
!      if( igas .ne. 0 )then
!        if( carma%f_element(ielem)%f_itype .eq. I_COREMASS .and. carma%f_element(ielem)%f_isolute .eq. 0 )then
!          if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - inucgas ne 0 but isolute eq 0'
!          rc = RC_ERROR
!          return
!        endif
!      endif
    enddo

    do ielem = 1, carma%f_NELEM
      if( carma%f_nnuc2elem(ielem) .gt. 0 ) then
        do inuc2 = 1, carma%f_nnuc2elem(ielem)
          if( carma%f_inuc2elem(inuc2,ielem) .gt. carma%f_NELEM )then
            if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - component of inuc2elem > NELEM'
            rc = RC_ERROR
            return
          endif
        enddo
      endif
    enddo

    ! Particle grids are incompatible if there is no target bin with enough
    ! mass to accomodate nucleated particle.
    bad_grid = .false.

    do iefrom = 1,carma%f_NELEM   ! source element

      igfrom = carma%f_element(iefrom)%f_igroup
      neto   = carma%f_nnuc2elem(iefrom)

      if( neto .gt. 0 )then

        do inuc2 = 1,neto
          ieto = carma%f_inuc2elem(inuc2,iefrom)
          igto = carma%f_element(ieto)%f_igroup

          do ifrom = 1,carma%f_NBIN   ! source bin
            if( carma%f_inuc2bin(ifrom,igfrom,igto) .eq. 0 )then
              if (carma%f_do_print) write(carma%f_LUNOPRT,7) igfrom,ifrom,igto
              bad_grid = .true.
            endif
          enddo
        enddo
      endif
    enddo

    if( bad_grid )then
      if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - incompatible grids for nucleation'
      rc = RC_ERROR
      return
    endif
      
    if (carma%f_do_print_init) then
    
      ! Report some initialization values!
      write(carma%f_LUNOPRT,5)
      write(carma%f_LUNOPRT,1) 'inucgas  ',(carma%f_inucgas(i),i=1,carma%f_NGROUP)
      write(carma%f_LUNOPRT,1) 'inuc2elem',(carma%f_inuc2elem(1,i),i=1,carma%f_NELEM)
      write(carma%f_LUNOPRT,1) 'ievp2elem',(carma%f_ievp2elem(i),i=1,carma%f_NELEM)
      write(carma%f_LUNOPRT,1) 'isolute ',(carma%f_element(i)%f_isolute,i=1,carma%f_NELEM)
    
      do isol = 1,carma%f_NSOLUTE
        write(carma%f_LUNOPRT,2) 'solute number   ',isol
        write(carma%f_LUNOPRT,3) 'solute name:    ',carma%f_solute(isol)%f_name
        write(carma%f_LUNOPRT,4) 'molecular weight',carma%f_solute(isol)%f_wtmol
        write(carma%f_LUNOPRT,4) 'mass density    ',carma%f_solute(isol)%f_rho    
      enddo
    endif
    
    
    ! Initialize indexes for the gases and check to make sure if H2SO4 is used
    ! that it occurs after H2O. This is necessary for supersaturation calculations.
    carma%f_igash2o   = 0
    carma%f_igash2so4 = 0
    carma%f_igasso2   = 0

    do igas = 1, carma%f_NGAS
      if (carma%f_gas(igas)%f_icomposition == I_GCOMP_H2O) then
        carma%f_igash2o = igas
      else if (carma%f_gas(igas)%f_icomposition == I_GCOMP_H2SO4) then
        carma%f_igash2so4 = igas
      else if (carma%f_gas(igas)%f_icomposition == I_GCOMP_SO2) then
        carma%f_igasso2 = igas
      end if
    end do
    
    if ((carma%f_igash2so4 /= 0) .and. (carma%f_igash2o > carma%f_igash2so4)) then
      if (carma%f_do_print) write(carma%f_LUNOPRT,*) 'CARMA_InitializeGrowth::ERROR - H2O gas must come before H2SO4.'
      rc = RC_ERROR
      return
    end if

    return
  end subroutine CARMA_InitializeGrowth         

  !! Calculate the optical properties for each particle bin at each of
  !! the specified wavelengths. The optical properties include the
  !! extinction efficiency, the single scattering albedo and the
  !! asymmetry factor.
  !!
  !! NOTE: For these calculations, the particles are assumed to be spheres and
  !! Mie code is used to calculate the optical properties.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeOptics(carma, rc)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
    
    integer                            :: igroup      ! group index
    integer                            :: iwave       ! wavelength index
    integer                            :: ibin        ! bin index
    real(kind=f)                       :: Qext
    real(kind=f)                       :: Qsca
    real(kind=f)                       :: asym
   
    
    ! Assume success.
    rc = RC_OK    
    
    ! Were any wavelengths specified?
    do iwave = 1, carma%f_NWAVE
      do igroup = 1, carma%f_NGROUP
     
        ! Should we calculate mie properties for this group?
        if (carma%f_group(igroup)%f_do_mie) then 
       
          do ibin = 1, carma%f_NBIN

            ! Assume the particle is homogeneous (no core).
            !
            ! NOTE: The miess does not converge over as broad a
            ! range of input parameters as bhmie, but it can handle
            ! coated spheres.

            call mie(carma, &
                     carma%f_group(igroup)%f_imiertn, &
                     carma%f_group(igroup)%f_r(ibin), &
                     carma%f_wave(iwave), &
                     carma%f_group(igroup)%f_refidx(iwave), &
                     Qext, &
                     Qsca, &
                     asym, &
                     rc)

            if (rc < RC_OK) then
              if (carma%f_do_print) then
                 write(carma%f_LUNOPRT, *) "CARMA_InitializeOptics:: &
                      &Mie failed for (band, wavelength, group, bin)", &
                      iwave, carma%f_wave(iwave), igroup, ibin
              end if
              return
            end if

            carma%f_group(igroup)%f_qext(iwave, ibin) = Qext
            carma%f_group(igroup)%f_ssa(iwave, ibin)  = Qsca / Qext
            carma%f_group(igroup)%f_asym(iwave, ibin) = asym

          end do
        end if
      end do
    end do   
    
    return
  end subroutine CARMA_InitializeOptics         

  !! Perform initialization of variables related to thermodynamical calculations that
  !! are not dependent on the model state.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeThermo(carma, rc)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
        
    ! Assume success.
    rc = RC_OK

    return
  end subroutine CARMA_InitializeThermo         

  !! Perform initialization of variables related to vertical transport that are not dependent
  !! on the model state.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeVertical(carma, rc, vf_const)
    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
    real(kind=f), intent(in), optional :: vf_const
        
    ! Assume success.
    rc = RC_OK

    ! Was a constant vertical velocity specified?
    carma%f_ifall = 1
    carma%f_vf_const = 0._f
    
    if (present(vf_const)) then
      if (vf_const /= 0._f) then
        carma%f_ifall = 0
        carma%f_vf_const = vf_const
      end if
    end if
    
    ! Specify the boundary conditions for vertical transport.
    carma%f_itbnd_pc  = I_FIXED_CONC
    carma%f_ibbnd_pc  = I_FIXED_CONC
    
    return
  end subroutine CARMA_InitializeVertical         

  !! The routine should be called when the carma object is no longer needed. It deallocates
  !! any memory allocations made by CARMA (during CARMA_Create()), and failure to call this
  !!routine could result in memory leaks.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_Create
  subroutine CARMA_Destroy(carma, rc)
    use carmaelement_mod
    use carmagas_mod
    use carmagroup_mod
    use carmasolute_mod

    type(carma_type), intent(inout)    :: carma
    integer, intent(out)               :: rc
    
    ! Local variables
    integer   :: ier
    integer   :: igroup
    integer   :: ielem
    integer   :: isolute
    integer   :: igas
    
    ! Assume success.
    rc = RC_OK
    
    ! If allocated, deallocate all the variables that were allocated in the Create() method.
    if (allocated(carma%f_group)) then
      do igroup = 1, carma%f_NGROUP
        call CARMAGROUP_Destroy(carma, igroup, rc)
        if (rc < 0) return
      end do
      
      deallocate( &
        carma%f_group, &
        carma%f_icoag, &
        carma%f_inucgas, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Destroy: ERROR deallocating groups, status=", ier
        rc = RC_ERROR
      endif
    endif
    
    if (allocated(carma%f_element)) then
      do ielem = 1, carma%f_NELEM
        call CARMAELEMENT_Destroy(carma, ielem, rc)
        if (rc < RC_OK) return
      end do
      
      deallocate( &
        carma%f_element, &
        carma%f_igrowgas, &
        carma%f_inuc2elem, &
        carma%f_inucproc, &
        carma%f_ievp2elem, &
        carma%f_nnuc2elem, &
        carma%f_nnucelem, &
        carma%f_inucelem, &
        carma%f_if_nuc, &
        carma%f_rlh_nuc, &
        carma%f_icoagelem, &
        carma%f_icoagelem_cm, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Destroy: ERROR deallocating elements, status=", ier
        rc = RC_ERROR
      endif
    endif
    
    if (allocated(carma%f_inuc2bin)) then
      deallocate( &
        carma%f_inuc2bin, &
        carma%f_ievp2bin, &
        carma%f_nnucbin, &
        carma%f_inucbin, &
        carma%f_diffmass, &
        carma%f_volx, &
        carma%f_ilow, &
        carma%f_jlow, &
        carma%f_iup, &
        carma%f_jup, &
        carma%f_npairl, &
        carma%f_npairu, &
        carma%f_iglow, &
        carma%f_jglow, &
        carma%f_igup, &
        carma%f_jgup, &
        carma%f_kbin, &
        carma%f_pkernel, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Destroy: ERROR deallocating bins, status=", ier
        rc = RC_ERROR
      endif
    endif
    
    if (carma%f_NSOLUTE > 0) then
      do isolute = 1, carma%f_NSOLUTE
        call CARMASOLUTE_Destroy(carma, isolute, rc)
        if (rc < RC_OK) return
      end do
      
      if (allocated(carma%f_solute)) then
        deallocate( &
          carma%f_solute, &
          stat=ier) 
        if(ier /= 0) then
          if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Destroy: ERROR deallocating solutes, status=", ier 
          rc = RC_ERROR
        endif
      endif
    end if
     
    if (carma%f_NGAS > 0) then
      do igas = 1, carma%f_NGAS
        call CARMAGAS_Destroy(carma, igas, rc)
        if (rc < RC_OK) return
      end do
      
      if (allocated(carma%f_gas)) then
        deallocate( &
          carma%f_gas, &
          stat=ier) 
        if(ier /= 0) then
          if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Destroy: ERROR deallocating gases, status=", ier
          rc = RC_ERROR
        endif
      endif
    end if     
    
    if (carma%f_NWAVE > 0) then
      if (allocated(carma%f_wave)) then
        deallocate( &
          carma%f_wave, &
          carma%f_dwave, &
          carma%f_do_wave_emit, &
          stat=ier) 
        if(ier /= 0) then
          if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_Destroy: ERROR deallocating wavelengths, status=", ier
          rc = RC_ERROR
          return
        endif
      endif
    endif
    
    return
  end subroutine CARMA_Destroy

  ! Configuration    
        
  !! Add a coagulation process between two groups (<i>igroup1</i> and <i>igroup2</i>), with the resulting
  !! particle being in the destination group (<i>igroup3</i>). If <i>ck0</i> is specifed, then a constant
  !! coagulation kernel will be used.
  subroutine CARMA_AddCoagulation(carma, igroup1, igroup2, igroup3, icollec, rc, ck0, grav_e_coll0)
    type(carma_type), intent(inout)    :: carma         !! the carma object
    integer, intent(in)                :: igroup1       !! first source group
    integer, intent(in)                :: igroup2       !! second source group
    integer, intent(in)                :: igroup3       !! destination group
    integer, intent(in)                :: icollec       !! collection technique [I_COLLEC_CONST | I_COLLEC_FUCHS | I_COLLEC_DATA] 
    integer, intent(out)               :: rc            !! return code, negative indicates failure
    real(kind=f), intent(in), optional :: ck0           !! if specified, forces a constant coagulation kernel
    real(kind=f), intent(in), optional :: grav_e_coll0  !! if <i>icollec</i> is I_COLLEC_CONST
                                                        !! the constant gravitational collection efficiency
    
    ! Assume success.
    rc = RC_OK

    ! Make sure the groups exists.
    if (igroup1 > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, '(a,i3,a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed group (", &
        igroup1, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if

    if (igroup2 > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, '(a,i3,a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed group (", &
        igroup2, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    if (igroup3 > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, '(a,i3,a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed group (", &
        igroup3, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Indicate that the groups coagulate together.
    carma%f_icoag(igroup1, igroup2) = igroup3
    
    ! If ck0 was specified, then we use a fixed coagulation rate of ck0.
    if (present(ck0)) then
      carma%f_ck0 = ck0
      carma%f_icoagop = I_COAGOP_CONST
    else
      carma%f_icoagop = I_COAGOP_CALC
    end if
    
    ! What collection technique is specified.
    if (icollec > I_COLLEC_DATA) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, '(a,i3,a)') "CARMA_AddCoagulation:: ERROR - The specifed collection method (", &
        icollec, ") is unknown."
      rc = RC_ERROR
      return
    end if
    
    if (icollec == I_COLLEC_CONST) then
      if (present(grav_e_coll0)) then
        carma%f_grav_e_coll0 = grav_e_coll0
      else
        if (carma%f_do_print) then
           write(carma%f_LUNOPRT, *) "CARMA_AddCoagulation::&
                &ERROR - A constant gravitational collection was requested, &
                &but grav_e_coll0 was not provided."
        end if
        rc = RC_ERROR
        return
      end if
    end if
    
    carma%f_icollec = icollec
    
    return
  end subroutine CARMA_AddCoagulation
    
  !! Add a growth process between the element (<i>ielem</i>) and gas (<i>igas</i>) specifed. The element
  !! and gas should have already been defined using <i>CARMA_AddElement()</i> and <i>CARMA_AddGas()</i>.
  !!
  !! NOTE: Each element can only have one volatile component.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_AddElement
  !! @see CARMA_AddGas
  subroutine CARMA_AddGrowth(carma, ielem, igas, rc)
    type(carma_type), intent(inout)    :: carma    !! the carma object
    integer, intent(in)                :: ielem    !! the element index
    integer, intent(in)                :: igas     !! the gas index
    integer, intent(out)               :: rc       !! return code, negative indicates failure
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure the element exists.
    if (ielem > carma%f_NELEM) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if

    ! Make sure there are enough gases allocated.
    if (igas > carma%f_NGAS) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%f_NGAS, ")."
      rc = RC_ERROR
      return
    end if
    
    ! If not already defined, indicate that the element can grow with the specified gas.
    if (carma%f_igrowgas(ielem) /= 0) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddGrowth:: ERROR - The specifed element (", &
        ielem, ") already has gas (", carma%f_igrowgas(ielem), ") condensing on it."
      rc = RC_ERROR
      return
    else 
      carma%f_igrowgas(ielem) = igas
    end if
    
    return
  end subroutine CARMA_AddGrowth
    
  !! Add a nucleation process that nucleates one element (<i>elemfrom</i>) to another element (<i>elemto</i>)
  !! using the specified gas (<i>igas</i>). The elements and gas should have already been defined
  !! using <i>CARMA_AddElement()</i> and <i>CARMA_AddGas()</i>. The nucleation scheme is indicated by
  !! inucproc, and can be one of:
  !!
  !!   - <i>I_DROPACT</i>
  !!   - <i>I_AERFREEZE</i>
  !!   - <i>I_DROPFREEZE</i>
  !!   - <i>I_ICEMELT</i>
  !!   - <i>I_HETNUC</i>
  !!   - <i>I_HOMNUC</i>
  !! 
  !! There are multiple parameterizations for I_AERFREEZE, so when that is selected the
  !! particular parameterization needs to be indicated by adding it to I_AERFREEZE. The
  !! specific routines are:
  !!
  !!   - <i>I_AF_TABAZADEH_2000</i>
  !!   - <i>I_AF_KOOP_2000</i>
  !!   - <i>I_AF_MOHLER_2010</i>
  !!   - <i>I_AF_MURRAY_2010</i>
  !!
  !! One or more of these routines may be selected, but in general one of the first
  !! three should be selected and then it can optionally be combined with the glassy
  !! aerosols (I_AF_MURRAY_2010).
  !! 
  !! Total evaporation transfers particle mass from the destination element back to the
  !! element indicated by <i>ievp2elem</i>. This relationship is not automatically generated,
  !! because multiple elements can nucleate to a particular element and therefore the 
  !! reverse mapping is not unique.
  !!
  !! NOTE: The gas used for nucleation must be the same for all nucleation defined from
  !! elements of the same group.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see I_DROPACT
  !! @see I_AERFREEZE
  !! @see I_DROPFREEZE
  !! @see I_ICEMELT
  !! @see I_HETNUC
  !! @see I_HOMNUC
  !! @see I_AF_TABAZADEH_2000
  !! @see I_AF_KOOP_2000
  !! @see I_AF_MOHLER_2010
  !! @see I_AF_MURRAY_2010
  !! @see CARMA_AddElement
  !! @see CARMA_AddGas
  subroutine CARMA_AddNucleation(carma, ielemfrom, ielemto, inucproc, &
      rlh_nuc, rc, igas, ievp2elem)
      
    use carmaelement_mod, only         : CARMAELEMENT_Get
    
    type(carma_type), intent(inout)    :: carma       !! the carma object
    integer, intent(in)                :: ielemfrom   !! the source element
    integer, intent(in)                :: ielemto     !! the destination element
    integer, intent(in)                :: inucproc    !! the nucleation process
                                                      !! [I_DROPACT | I_AERFREEZE | I_ICEMELT | I_HETNUC | I_HOMNUC]
    real(kind=f), intent(in)           :: rlh_nuc     !! the latent heat of nucleation [cm<sup>2</sup>/s<sup>2</sup>]
    integer, intent(out)               :: rc          !! return code, negative indicated failure
    integer, optional, intent(in)      :: igas        !! the gas
    integer, optional, intent(in)      :: ievp2elem   !! the element created upon evaporation
    
    integer                            :: igroup      !! group for source element
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure the elements exist.
    if (ielemfrom > carma%f_NELEM) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed element (", &
        ielemfrom, ") is larger than the number of elements (", carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if

    if (ielemto > carma%f_NELEM) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed element (", &
        ielemto, ") is larger than the number of elements (", carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if

    if (present(ievp2elem)) then
      if (ievp2elem > carma%f_NELEM) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed element (", &
          ievp2elem, ") is larger than the number of elements (", carma%f_NELEM, ")."
        rc = RC_ERROR
        return
      end if
    end if


    ! Make sure there are enough gases allocated.
    if (present(igas)) then
      if (igas > carma%f_NGAS) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMA_AddNucleation:: ERROR - The specifed gas (", &
          igas, ") is larger than the number of gases (", carma%f_NGAS, ")."
        rc = RC_ERROR
        return
      end if
    end if
    
    
    ! If aerosol freezing is selected, but no I_AF_xxx sub-method is selected, then indicate an error.
    if (inucproc == I_AERFREEZE) then
      if (carma%f_do_print) then
         write(carma%f_LUNOPRT, *) "CARMA_AddNucleation::&
              &ERROR - I_AERFREEZE was specified without an I_AF_xxx value."
      end if
      return
    end if
    
    
    ! Array <inucgas> maps a particle group to its associated gas for nucleation:
    ! Nucleation from group <igroup> is associated with gas <inucgas(igroup)>
    ! Set to zero if particles are not subject to nucleation.
    if (present(igas)) then
      call CARMAELEMENT_Get(carma, ielemfrom, rc, igroup=igroup)
      
      if (rc >= RC_OK) then
        carma%f_inucgas(igroup) = igas
      end if
    end if


    ! Nucleation transfers particle mass from element <ielem> to element
    ! <inuc2elem(i,ielem)>, where <i> ranges from 0 to the number of elements
    !  nucleating from <ielem>.
!    carma%f_nnucelem(ielemto) = carma%f_nnucelem(ielemto) + 1
!    carma%f_inucelem(carma%f_nnucelem(ielemto), ielemto) = ielemfrom
    carma%f_nnuc2elem(ielemfrom) = carma%f_nnuc2elem(ielemfrom) + 1
    carma%f_inuc2elem(carma%f_nnuc2elem(ielemfrom), ielemfrom) = ielemto
!    carma%f_if_nuc(ielemfrom,carma%f_inuc2elem(carma%f_nnuc2elem(ielemfrom), ielemfrom)) = .true.

    ! <inucproc(iefrom,ieto)> specifies what nucleation process nucleates
    ! particles from element <ielem> to element <ieto>:
    !   I_DROPACT:  Aerosol activation to droplets 
    !   I_AERFREEZE: Aerosol homogeneous freezing
    !   I_DROPFREEZE: Droplet homogeneous freezing
    !   I_GLFREEZE: Glassy Aerosol heteroogeneous freezing
    !   I_GLAERFREEZE: Glassy & Aerosol freezing
    carma%f_inucproc(ielemfrom, ielemto) = inucproc


    ! Total evaporation mapping: total evaporation transfers particle mass from
    ! element <ielem> to element <ievp2elem(ielem)>.
    !
    ! NOTE: This array is not automatically derived from <inuc2elem> because multiple
    ! elements can nucleate to a particular element (reverse mapping is not
    ! unique).
    if (present(ievp2elem)) carma%f_ievp2elem(ielemto) = ievp2elem


    ! <rlh_nuc(iefrom,ieto)> is the latent heat released by nucleation
    ! from element <iefrom> to element <ieto> [cm^2/s^2].
    carma%f_rlh_nuc(ielemfrom,ielemto) = rlh_nuc

    return
  end subroutine


  ! Query, Control and State I/O

  !! Gets the information about the carma object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_Create
  subroutine CARMA_Get(carma, rc, LUNOPRT, NBIN, NELEM, NGAS, NGROUP, NSOLUTE, NWAVE, do_detrain, &
    do_drydep, do_fixedinit, do_grow, do_print, do_print_init, do_thermo, wave, dwave, do_wave_emit)
    
    type(carma_type), intent(in)        :: carma                !! the carma object
    integer, intent(out)                :: rc                   !! return code, negative indicates failure
    integer, optional, intent(out)      :: NBIN                 !! number of radius bins per group
    integer, optional, intent(out)      :: NELEM                !! total number of elements
    integer, optional, intent(out)      :: NGROUP               !! total number of groups
    integer, optional, intent(out)      :: NSOLUTE              !! total number of solutes
    integer, optional, intent(out)      :: NGAS                 !! total number of gases
    integer, optional, intent(out)      :: NWAVE                !! number of wavelengths
    integer, optional, intent(out)      :: LUNOPRT              !! logical unit number for output
    logical, optional, intent(out)      :: do_detrain           !! do detrainement?
    logical, optional, intent(out)      :: do_drydep            !! do dry deposition?
    logical, optional, intent(out)      :: do_fixedinit         !! do initialization from reference atm?
    logical, optional, intent(out)      :: do_grow              !! do condensational growth?
    logical, optional, intent(out)      :: do_print             !! do print output?
    logical, optional, intent(out)      :: do_print_init        !! do print initialization output?
    logical, optional, intent(out)      :: do_thermo            !! do thermodynamics?
    real(kind=f), optional, intent(out) :: wave(carma%f_NWAVE)  !! the wavelengths centers (cm)
    real(kind=f), optional, intent(out) :: dwave(carma%f_NWAVE) !! the wavelengths widths (cm)
    logical, optional, intent(out)      :: do_wave_emit(carma%f_NWAVE) !! do emission in this band?
    
    ! Assume success.
    rc = RC_OK
    
    if (present(LUNOPRT))  LUNOPRT = carma%f_LUNOPRT
    if (present(NBIN))     NBIN    = carma%f_NBIN
    if (present(NELEM))    NELEM   = carma%f_NELEM
    if (present(NGAS))     NGAS    = carma%f_NGAS
    if (present(NGROUP))   NGROUP  = carma%f_NGROUP
    if (present(NSOLUTE))  NSOLUTE = carma%f_NSOLUTE
    if (present(NWAVE))    NWAVE   = carma%f_NWAVE
    
    if (present(do_detrain))    do_detrain     = carma%f_do_detrain
    if (present(do_drydep))     do_drydep      = carma%f_do_drydep
    if (present(do_grow))       do_grow        = carma%f_do_grow
    if (present(do_fixedinit))  do_fixedinit   = carma%f_do_fixedinit
    if (present(do_print))      do_print       = carma%f_do_print
    if (present(do_print_init)) do_print_init  = carma%f_do_print_init
    if (present(do_thermo))     do_thermo      = carma%f_do_thermo

    if (present(wave))  wave(:)    = carma%f_wave(:)
    if (present(dwave)) dwave(:)   = carma%f_dwave(:)
    if (present(do_wave_emit)) do_wave_emit(:)   = carma%f_do_wave_emit(:)

    return
  end subroutine CARMA_Get

end module
