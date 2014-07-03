!! This module defines types used in the CARMA module. The types need to be defined here
!! to avoid circular references between different modules (e.g. carma_mod and
!! carmastate_mod).
!!
!! NOTE: All the field members are prefixed by f_. This is done because of the macros that
!! are used to map between the older F77 common block names for variables to the newer F90
!! structure member names for the fields. This is done in carma_globaer.h to keep the core
!! CARMA code looking similar to the F77 code to make it easier for scientists with CARMA
!! experience to port their code.  Some compilers (e.g Portland Group) have preprocessors
!! that will fail to handle the macros in carma_globaer.h properly resulting in recursion
!! errors during compiling. By making the field member name different, the recursion
!! problems should be avoided.
!!
!! @version July-2009 
!! @author  Chuck Bardeen 
module carma_types_mod
  use carma_precision_mod
  use carma_constants_mod

  !! The CARMAELEMENT data type represents one of the components of a cloud or aerosol particle.
  !!
  !! The procedure for adding a variable to the CARMAELEMENT data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmaelement_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMAELEMENT_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMAELEMENT data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmaelement_type
  
    !   name          Name of the element
    !   shortname     Short name of the element
    !   rho           Mass density of particle element [g/cm^3]
    !   igroup        Group to which the element belongs
    !   itype         Particle type specification
    !   icomposition  Particle compound specification
    !   isolute       Index of solute for the particle element
    !
    character(len=CARMA_NAME_LEN)               :: f_name
    character(len=CARMA_SHORT_NAME_LEN)         :: f_shortname
    real(kind=f), allocatable, dimension(:)     :: f_rho          ! (NBIN)
    integer                                     :: f_igroup
    integer                                     :: f_itype
    integer                                     :: f_icomposition
    integer                                     :: f_isolute
  end type carmaelement_type


  !! The CARMAGAS data type represents a gas.
  !!
  !! The procedure for adding a variable to the CARMAGAS data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmagas_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMAGAS_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMAGAS data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmagas_type
  
    !   name          Name of the gas
    !   shortname     Short name of the gas
    !   wtmol         Molecular weight for the gas [g/mol]
    !   ivaprtn       vapor pressure routine for the gas
    !   dgc_threshold convergence criteria for gas concentration [fraction]
    !   ds_threshold  convergence criteria for gas saturation [fraction]
    !
    character(len=CARMA_NAME_LEN)               :: f_name
    character(len=CARMA_SHORT_NAME_LEN)         :: f_shortname
    real(kind=f)                                :: f_wtmol
    integer                                     :: f_ivaprtn
    integer                                     :: f_icomposition
    real(kind=f)                                :: f_dgc_threshold
    real(kind=f)                                :: f_ds_threshold
  end type carmagas_type


  !! The CARMAGROUP data type represents a cloud or aerosol partcile.
  !!
  !! The procedure for adding a variable to the CARMAGROUP data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmagroup_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMAGROUP_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMAGROUP data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmagroup_type
  
    !   name        Name of the particle
    !   shortname   Short name of the particle
    !   cnsttype    constituent type [I_CNSTTYPE_PROGNOSTIC | I_CNSTTYPE_DIAGNOSTIC]
    !   maxbin      the last prognostic bin in the group
    !   nelem       Number of elements in group           
    !   ncore       Number of core elements (itype = 2) in group           
    !   ishape      Describes particle shape for group
    !   ienconc     Particle number conc. element for group
    !   imomelem    Scondary moment element for group
    !   icorelem    Core elements (itype = 2) in group           
    !   solfac      Solubility factor for wet deposition
    !   is_ice      If .true. then ice particle 
    !   is_cloud    If .true. then cloud particle
    !   is_sulfate  If .true. then sulfate particle
    !   do_mie      If .true. then do mie calculations 
    !   do_wetdep   If .true. then do wet deposition 
    !   grp_do_drydep If .true. then do dry deposition 
    !   grp_do_vtran If .true. then do sedimentation 
    !   scavcoef    Scavenging coefficient for wet deopistion (1/mm)
    !   if_sec_mom  If .true. then core second moment (itype = 3) used    {setupgrow}
    !   irhswell    Indicates method for swelling particles from RH
    !   irhswcomp   Indicates composition for swelling particles from RH
    !   rmin        Radius of particle in first bin [cm]
    !   rmassmin    Mass of particle in first bin [g]
    !   rmrat       Ratio of masses of particles in consecutive bins
    !   eshape      Ratio of particle length / diameter 
    !   r           Radius bins [cm]
    !   rmass       Mass bins [g]
    !   rrat        Ratio of maximum diameter to diameter of equivalent sphere
    !   arat        Ratio of projected area to projected area of containing sphere
    !   vol         Particle volume [cm^3]
    !   dr          Width of bins in radius space [cm]
    !   dm          Width of bins in mass space [g]
    !   rmassup     Upper bin boundary mass [g]
    !   rup         Upper bin boundary radius [cm]
    !   rlow        Lower bin boundary radius [cm]
    !   refidx      refractive index
    !   qext        extinction efficiency
    !   ssa         single scattering albedo
    !   asym        asymmetry factor
    !   ifallrtn    routine to use to calculate fall velocity  [I_FALLRTN_...]
    !   imiertn     mie routine for optical properties [I_MIERTN_...]
    !   dpc_threshold convergence criteria for particle concentration [fraction]
    !
    character(len=CARMA_NAME_LEN)               :: f_name
    character(len=CARMA_SHORT_NAME_LEN)         :: f_shortname
    integer                                     :: f_cnsttype
    integer                                     :: f_maxbin
    integer                                     :: f_nelem
    integer                                     :: f_ncore
    integer                                     :: f_ishape
    integer                                     :: f_ienconc
    integer                                     :: f_imomelem
    real(kind=f)                                :: f_solfac
    real(kind=f)                                :: f_scavcoef
    logical                                     :: f_if_sec_mom
    logical                                     :: f_is_ice
    logical                                     :: f_is_cloud
    logical                                     :: f_is_sulfate
    logical                                     :: f_do_mie
    logical                                     :: f_do_wetdep
    logical                                     :: f_grp_do_drydep
    logical                                     :: f_grp_do_vtran
    integer                                     :: f_irhswell
    integer                                     :: f_irhswcomp
    integer                                     :: f_ifallrtn
    integer                                     :: f_imiertn
    real(kind=f)                                :: f_rmin
    real(kind=f)                                :: f_rmassmin
    real(kind=f)                                :: f_rmrat
    real(kind=f)                                :: f_eshape
    real(kind=f), allocatable, dimension(:)     :: f_r          ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_rmass      ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_vol        ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_dr         ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_dm         ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_rmassup    ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_rup        ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_rlow       ! (NBIN)
    complex(kind=f), allocatable, dimension(:)  :: f_refidx     ! (NWAVE)
    real(kind=f), allocatable, dimension(:,:)   :: f_qext       ! (NWAVE,NBIN)
    real(kind=f), allocatable, dimension(:,:)   :: f_ssa        ! (NWAVE,NBIN)
    real(kind=f), allocatable, dimension(:,:)   :: f_asym       ! (NWAVE,NBIN)
    integer, allocatable, dimension(:)          :: f_icorelem   ! (NELEM)
    real(kind=f), allocatable, dimension(:)     :: f_arat       ! (NBIN)
    real(kind=f), allocatable, dimension(:)     :: f_rrat       ! (NBIN)
    real(kind=f)                                :: f_dpc_threshold
  end type carmagroup_type
  
  
  !! The CARMASOLUTE data type represents a gas.
  !!
  !! The procedure for adding a variable to the CARMASOLUTE data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmagas_type is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMASOLUTE_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMASOLUTE data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! routines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carmasolute_type
  
    !   name        Name of the solute
    !   shortname   Short name of the solute
    !   ions        Number of ions solute dissociates into
    !   wtmol       Molecular weight of solute
    !   rho         Mass density of solute
    !
    character(len=CARMA_NAME_LEN)               :: f_name
    character(len=CARMA_SHORT_NAME_LEN)         :: f_shortname
    integer                                     :: f_ions
    real(kind=f)                                :: f_wtmol
    real(kind=f)                                :: f_rho
  end type carmasolute_type


  !! The CARMA data type replaces the common blocks that were used in the F77 version of
  !! CARMA. This allows the code to be written to allow for multiple threads to call CARMA
  !! routines simulataneously. This thread safety is necessary for to run CARMA under OPEN/MP.
  !!
  !! The procedure for adding a variable to the CARMA data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the appropriate create or initialization routine.
  !!    - Deallocate the variable in the approprate finalize and destroy routines.
  !!  - Add an alias for the variable to carma_globaer.h and associate it with the variable
  !!     in this typedef.  
  !!
  !! NOTE: While the carmatype is public, routines outside of the CARMA module should not look
  !! at or manuipulate fields of this structure directly. There should be CARMA_XXX methods
  !! to do anything that is needed with this structure, and use of these methods will allow
  !! the CARMA data type structure to evolve without impacting code in the parent model.
  !! The contents of the structure had to be made public, since the CARMA microphysics
  !! rountines are implemented in separate files outside of this model; however, logically
  !! they are part of the model and are the only routines outside of this module that should
  !! access fields of this structure directly.
  type, public :: carma_type
  
    ! Model Dimensions
    !
    !  NGROUP   number of particle groups
    !  NELEM    number of particle components (elements)
    !  NBIN     number of size bins per element
    !  NGAS     number of gases (may be 0)
    !  NSOLUTE  number of solutes (may be 0)
    !  NWAVE    number of wavelength bands (may be 0)
    !
    integer :: f_NGROUP
    integer :: f_NELEM
    integer :: f_NBIN
    integer :: f_NGAS
    integer :: f_NSOLUTE
    integer :: f_NWAVE

    ! Output logical unit numbers
    !
    ! NOTE: CARMA will not directly access files or keep track of file names. It is the
    ! parent model's responsibility to provide the logical unit number to be used for
    ! model output.
    !
    integer :: f_LUNOPRT   ! output print file

    ! Model startup control variables
    !
    !   do_print  .t. if print output is desired
    !
    logical            :: f_do_print
    
    
    ! Configuration Objects
    !
    ! These are all other objects that are parts of the CARMA model. This is
    ! an attempt to break up the large common block that has historically been
    ! the structure of CARMA so the code is easier to understand and to
    ! maintain.
    !
    !   element     Particle component    
    !   gas         Gas   
    !   group       Particle    
    !   solute      Element solute    
    !
    ! NOTE: In the future, it may make sense to create objects that represent
    ! the CARMA processes. This would encapsulate all the variables related to
    ! a particular process into one structure. Candidate processes include:
    ! transport, growth, nucleation, coagulation, ...
    !
    type(carmaelement_type), allocatable, dimension(:)    :: f_element    ! (NELEM)
    type(carmagas_type),     allocatable, dimension(:)    :: f_gas        ! (NGAS)
    type(carmagroup_type),   allocatable, dimension(:)    :: f_group      ! (NGROUP)
    type(carmasolute_type),  allocatable, dimension(:)    :: f_solute     ! (NSOLUTE)



    ! Model option & control variables
    !
    !   conmax      Minumum relative concentration to consider in varstep   {prestep}   
    !   icoag       Coagulation mapping array                           {setupcoag}
    !   icoagelem   Coagulation element mapping array                   {setupcoag}
    !   icoagelem_cm Coagulation element mapping array for second mom   {setupcoag}
    !   ifall       Fall velocity options                               {setupvfall}
    !   icoagop     Coagulation kernel options                          {setupckern}
    !   icollec     Gravitational collection options                      {setupckern}
    !   itbnd_pc    Top boundary condition flag for particles             {init}
    !   ibbnd_pc    Bottom boundary condition flag for particles          {init}
    !   do_vdiff    If .true. then do Brownian diffusion                  {init}
    !   do_coag     If .true. then do coagulation                         {init}
    !   do_detrain  If .true. then do detrainment                         {init}
    !   do_drydep   If .true. then do dry deposition                      {init}
    !   do_fixedinitIf .true. then do initialize from reference atm       {init}
    !   do_grow     If .true. then do condensational growth and evap.     {init}
    !   do_clearsky If .true. then do clear sky growth and coagulation    {init}
    !   do_incloud  If .true. then do incloud growth and coagulation      {init}
    !   do_explised If .true. then do sedimentation with substepping      {init}
    !   do_pheat    If .true. then do particle heating for growth rates   {init}
    !   do_pheatatm If .true. then do particle heating on atmosphere      {init}
    !   do_print_init If .true. then do print initializtion info          {init}
    !   do_step     if .true. then varstepping succeeded                  {init}
    !   do_substep  if .true. then use substepping                        {init}
    !   do_thermo   if .true. then do solve thermodynamic equation        {init}
    !   do_vdiff    If .true. then do Brownian diffusion                  {init}
    !   do_vtran    If .true. then do vertical transport                  {init}
    !   do_cnst_rlh If .true. then uses constants for rlhe and rlhm       {setupgrow}
    !   igrowgas    Gas that condenses into a particle element            {setupgrow}
    !   inucgas     Gas that nucleates a particle group                   {setupnuc}
    !   if_nuc      Nucleation conditional array                          {setupaer}
    !   inucproc    Nucleation conditional array                          {setupaer}
    !   nnuc2elem   Number of elements that nucleate to element           {setupnuc}
    !   inuc2elem   Nucleation transfers particles into element inuc2elem {setupnuc}
    !   ievp2elem   Total evap. transfers particles into group ievp2elem  {setupnuc}
    !   ievp2bin    Total evap. transfers particles into bin ievp2bin     {setupnuc}
    !   inuc2bin    Nucleation transfers particles into bin inuc2bin      {setupnuc}
    !   maxsubsteps Maximum number of time substeps allowed
    !   minsubsteps Maximum number of time substeps allowed
    !   maxretries  Maximum number of substepping retries allowed
    !   igash2o     gas index for H2O
    !   igash2so4   gas index for H2SO4
    !   igasso2     gas index for SO2
    !   dt_threshold convergence criteria for temperature [fraction]
    !   cstick      accommodation coefficient - coagulation
    !   gsticki     accommodation coefficient - growth (ice), default = 0.93
    !   gstickl     accommodation coefficient - growth (liquid), default = 1.0
    !   tstick      accommodation coefficient - temperature, default = 1.0
    !
    logical                                       :: f_do_vdiff
    logical                                       :: f_do_drydep
    logical                                       :: f_do_coag
    logical                                       :: f_do_detrain
    logical                                       :: f_do_fixedinit
    logical                                       :: f_do_grow
    logical                                       :: f_do_clearsky
    logical                                       :: f_do_incloud
    logical                                       :: f_do_vtran
    logical                                       :: f_do_explised
    logical                                       :: f_do_pheat
    logical                                       :: f_do_pheatatm
    logical                                       :: f_do_print_init
    logical                                       :: f_do_step
    logical                                       :: f_do_substep
    logical                                       :: f_do_thermo
    logical                                       :: f_do_cnst_rlh
    logical, allocatable, dimension(:,:)          :: f_if_nuc       !(NELEM,NELEM)
    real(kind=f)                                  :: f_conmax
    integer                                       :: f_igash2o
    integer                                       :: f_igash2so4
    integer                                       :: f_igasso2
    integer                                       :: f_maxsubsteps 
    integer                                       :: f_minsubsteps 
    integer                                       :: f_maxretries 
    integer                                       :: f_ifall
    integer                                       :: f_icoagop
    integer                                       :: f_icollec
    integer                                       :: f_itbnd_pc
    integer                                       :: f_ibbnd_pc
    integer, allocatable, dimension(:)            :: f_inucgas      ! NGROUP
    integer, allocatable, dimension(:)            :: f_igrowgas     ! NELEM
    integer, allocatable, dimension(:)            :: f_nnuc2elem    ! NELEM
    integer, allocatable, dimension(:)            :: f_ievp2elem    ! NELEM
    integer, allocatable, dimension(:)            :: f_nnucelem     ! NELEM
    integer, allocatable, dimension(:,:)          :: f_icoag        ! (NGROUP,NGROUP)
    integer, allocatable, dimension(:,:)          :: f_inucproc     ! (NELEM,NELEM)
    integer, allocatable, dimension(:,:)          :: f_inuc2elem    ! (NELEM,NELEM)
    integer, allocatable, dimension(:,:)          :: f_icoagelem    ! (NELEM,NGROUP)
    integer, allocatable, dimension(:,:)          :: f_icoagelem_cm ! (NELEM,NGROUP)
    integer, allocatable, dimension(:,:)          :: f_inucelem     ! (NELEM,NELEM*NGROUP)
    integer, allocatable, dimension(:,:,:)        :: f_inuc2bin     ! (NBIN,NGROUP,NGROUP)
    integer, allocatable, dimension(:,:,:)        :: f_ievp2bin     ! (NBIN,NGROUP,NGROUP)
    integer, allocatable, dimension(:,:,:)        :: f_nnucbin      ! (NGROUP,NBIN,NGROUP)
    integer, allocatable, dimension(:,:,:,:)      :: f_inucbin      ! (NBIN*NGROUP,NGROUP,NBIN,NGROUP)
    real(kind=f)                                  :: f_dt_threshold
    real(kind=f)                                  :: f_tstick
    real(kind=f)                                  :: f_gsticki
    real(kind=f)                                  :: f_gstickl
    real(kind=f)                                  :: f_cstick
  

    ! Particle bin structure
    !
    !   diffmass  Difference between <rmass> values
    !
    real(kind=f), allocatable, dimension(:,:,:,:)   :: f_diffmass   ! (NBIN,NGROUP,NBIN,NGROUP)

    !  Coagulation kernels and bin pair mapping
    !
    !   ck0           Constant coagulation kernel           {setupaer}
    !   grav_e_coll0  Constant value for collection effic.  {setupaer}
    !   volx          Coagulation subdivision variable      {setupcoag}
    !   ilow          Bin pairs for coagulation production  {setupcoag}
    !   jlow          Bin pairs for coagulation production  {setupcoag}
    !   iup           Bin pairs for coagulation production  {setupcoag}
    !   jup           Bin pairs for coagulation production  {setupcoag}
    !   npairl        Bin pair indices                      {setupcoag}
    !   npairu        Bin pair indices                      {setupcoag}
    !   kbin          lower bin for coagulation             {setupcoag}
    !   pkernel       Coagulation production variables      {setupcoag}
    !
    real(kind=f)                                        :: f_ck0
    real(kind=f)                                        :: f_grav_e_coll0
    real(kind=f), allocatable, dimension(:,:,:,:,:)     :: f_volx    ! (NGROUP,NGROUP,NGROUP,NBIN,NBIN)
    integer, allocatable, dimension(:,:,:)              :: f_ilow    ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:)              :: f_jlow    ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:)              :: f_iup     ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:)              :: f_jup     ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:)                :: f_npairl  ! (NGROUP,NBIN)
    integer, allocatable, dimension(:,:)                :: f_npairu  ! (NGROUP,NBIN)
    integer, allocatable, dimension(:,:,:,:,:)          :: f_kbin    ! (NGROUP,NGROUP,NGROUP,NBIN,NBIN)
    real(kind=f), allocatable, dimension(:,:,:,:,:,:)   :: f_pkernel ! (NBIN,NBIN,NGROUP,NGROUP,NGROUP,6)

    !  Coagulation group pair mapping
    !
    !   iglow      Group pairs for coagulation production  {setupcoag}
    !   jglow      Group pairs for coagulation production  {setupcoag}
    !   igup       Group pairs for coagulation production  {setupcoag}
    !   jgup       Group pairs for coagulation production  {setupcoag}
    !
    integer, allocatable, dimension(:,:,:) :: f_iglow  ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:) :: f_jglow  ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:) :: f_igup   ! (NGROUP,NBIN,NBIN*NBIN)
    integer, allocatable, dimension(:,:,:) :: f_jgup   ! (NGROUP,NBIN,NBIN*NBIN)

    !  Particle fall velocities
    !
    !   vf_const  Constant vertical fall velocity when ifall=0          {setupaer}
    !
    real(kind=f)                                        :: f_vf_const


    ! Condensational growth parameters
    !
    ! NOTE: Some of these variables are used for storing intermediate values in
    ! the calculations. They may no longer be necessary, when the code is
    ! implemented as F90 and values as passed as parameters between subroutines.
    !
    !   rlh_nuc   Latent heat released by nucleation [cm^2/s^2]       {setupaer}
    !   pratt     Terms in PPM advection scheme for condensation      {setupgkern}
    !   prat
    !   pden1
    !   palr
    real(kind=f), allocatable, dimension(:,:)        :: f_rlh_nuc    ! (NELEM,NELEM)
    real(kind=f), allocatable, dimension(:,:,:)      :: f_pratt      ! (3,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)      :: f_prat       ! (4,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)        :: f_pden1      ! (NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)        :: f_palr       ! (4,NGROUP)
    
    ! Optical Properties
    !   wave      Bin-center wavelengths [cm]
    !   dwave     width of radiation bands [cm]
    !   do_wave_emit If true, emission should be calculated the band
    !
    real(kind=f), allocatable, dimension(:)          :: f_wave       ! (NWAVE)
    real(kind=f), allocatable, dimension(:)          :: f_dwave      ! (NWAVE)
    logical, allocatable, dimension(:)               :: f_do_wave_emit  ! (NWAVE)
  end type carma_type  
 
 
  !! The cstate data type replaces portions of the common blocks that were used
  !! in the F77 version of CARMA. This allows the code to be written to allow for
  !! multiple threads to call CARMA routines simulataneously. This thread safety is
  !! necessary for to run CARMA under OPEN/MP.
  !!
  !! The procedure for adding a variable to the cstate data type is:
  !!  - Add the variable as a scalar or allocatable in the type definition.
  !!  - If the new variable is dynamic,
  !!    - Allocate the variable in the create routine.
  !!    - Deallocate the variable in the destroy routines.
  !!  - Add an alias for the variable to cstate.h and associate it with the
  !!      variable in this typedef.  
  !!
  !! NOTE: While the carmastate_type is public, routines outside of the CARMA module
  !! should not look at or manuipulate fields of this structure directly. There should
  !! be CARMASTATE_XXX methods to do anything that is needed with this structure, and
  !! use of these methods will allow the cstate data type structure to evolve without
  !! impacting code in the parent model. The contents of the structure had to be made
  !! public, since the CARMA microphysics rountines are implemented in separate files
  !! outside of this model; however, logically they are part of the model and are the
  !! only routines outside of this module that should access fields of this structure
  !! directly.
  type, public :: carmastate_type
  
    ! Parent CARMA object
    type(carma_type), pointer                   :: f_carma

    ! Model Dimensions
    !
    !  NZ       number of grid points in the column
    !  NZP1     NZ+1
    !  NGROUP   number of particle groups
    !  NELEM    number of particle components (elements)
    !  NBIN     number of size bins per element
    !  NGAS     number of gases (may be 0)
    !
    integer :: f_NZ
    integer :: f_NZP1
    
    ! Model option & control variables
    !
    !   time        Simulation time at end of current timestep [s]
    !   dtime       Substep Timestep size [s]
    !   dtime_orig  Original Timestep size [s]
    !   nretries    Number of substepping retries attempted
    real(kind=f)                                  :: f_time
    real(kind=f)                                  :: f_dtime
    real(kind=f)                                  :: f_dtime_orig
    real(kind=f)                                  :: f_nretries

    !   max_nretry  Maximum number of retries in a step
    !   nstep       Total number of steps taken
    !   nsubstep    Total number of substeps taken
    !   nretry      Total number of retries taken
    integer                                       :: f_max_nsubstep
    real(kind=f)                                  :: f_max_nretry
    real(kind=f)                                  :: f_nstep
    integer                                       :: f_nsubstep
    real(kind=f)                                  :: f_nretry
    
    real(kind=f), allocatable, dimension(:)       :: f_zsubsteps   ! (NZ)


    ! Model Grid
    !
    !  igridv     flag to specify desired vertical grid coord system    {initatm}
    !  igridh     flag to specify desired horizontal grid coord system  {initatm}
    !  xmet       Horizontal ds/dx (ds is metric distance)              {initatm}
    !  ymet       Horizontal ds/dy (ds is metric distance)              {initatm}
    !  zmet       Vertical ds/dz (ds is metric distance)                {initatm}
    !  zmetl      Vertical ds/dz at edges (ds is metric distance)       {initatm}
    !  xc         Horizontal position at center of box                  {initatm}
    !  yc         Horizontal position at center of box                  {initatm}
    !  zc         Altitude at layer mid-point                           {initatm}
    !  dx         Horizontal grid spacing                               {initatm}
    !  dy         Horizontal grid spacing                               {initatm}
    !  dz         Thickness of vertical layers                          {initatm}
    !  zl         Altitude at top of layer                              {initatm}
    !  lon        Longitude [deg] at xc, yc                             {initatm}
    !  lat        Latitude [deg] at xc, yc                              {initatm}
    !
    integer :: f_igridv
    integer :: f_igridh
    real(kind=f), allocatable, dimension(:)     :: f_xmet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_ymet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_zmet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_zmetl  ! (NZP1)
    real(kind=f), allocatable, dimension(:)     :: f_xc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_yc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_zc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_dx     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_dy     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_dz     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_zl     ! (NZP1)
    real(kind=f)                                :: f_lon
    real(kind=f)                                :: f_lat

    ! Particle bin structure
    !
    !   rhop      Mass density of particle groups [g/cm^3]
    !   r_wet     Wet particle radius from RH swelling [cm]             {setupvfall}
    !   rlow_wet  Wet particle radius (lower bound) from RH swelling [cm]             {setupvfall}
    !   rup_wet   Wet particle radius (upper bound) from RH swelling [cm]             {setupvfall}
    !   rhop_wet  Wet Mass density of particle groups [g/cm^3]
    !   r_ref     Reference wet particle radius from RH swelling [cm]             {setupvfall}
    !   rhop_ref  Reference wet Mass density of particle groups [g/cm^3]
    !
    real(kind=f), allocatable, dimension(:,:,:) :: f_rhop       ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: f_rhop_wet   ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: f_r_wet      ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: f_rlow_wet   ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: f_rup_wet    ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: f_r_ref      ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:) :: f_rhop_ref   ! (NZ,NBIN,NGROUP)

    ! Primary model state variables
    !
    !  pc          Particle concentration [/x_units/y_units/z_units]  {initaer}
    !  pcd         Detrained particle concentration [/x_units/y_units/z_units]  {initaer}
    !  pc_surf     Particles on surface [/cm2]                        {initaer}
    !  sedimentationflux     Particles sedimented to surface [/cm2/s]                        {initaer}
    !  gc          Gas concentration [g/x_units/y_units/z_units]      {initgas}
    !  cldfrc      Cloud fraction [fraction]
    !  rhcrit      Relative humidity for onset of liquid clouds [fraction]
    !
    real(kind=f), allocatable, dimension(:,:,:) :: f_pc         ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:,:) :: f_pcd        ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)   :: f_pc_surf    ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)   :: f_sedimentationflux    ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)   :: f_gc         ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:)     :: f_cldfrc     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_rhcrit     ! (NZ)

    ! Secondary model variables
    !
    ! NOTE: Some of these variables are used for storing intermediate values in
    ! the calculations. They may no longer be necessary, when the code is
    ! implemented as F90 and values as passed as parameters between subroutines.
    !
    !   pcl         Particle concentration at beginning of time-step
    !   pconmax     Maximum particle concentration for each grid point
    !   gcl         Gas concentration at beginning of time-step
    !   d_gc        Change in gas concentration due to transport
    !   d_t        Change in temperature due to transport
    !   dpc_sed     Change in particle concentration due to sedimentation
    !   coaglg      Total particle loss rate due to coagulation for group
    !   coagpe      Particle production due to coagulation 
    !   rnuclg      Total particle loss rate due to nucleation for group
    !   rnucpe      Particle production due to nucleation 
    !   rhompe      Particle production due to homogeneous nucleation 
    !   pc_nucl     Particles produced due to nucleation (for the whole step, not just the substep)
    !   growlg      Total particle loss rate due to growth for group
    !   growle      Partial particle loss rate due to growth for element 
    !   growpe      Particle production due to growth 
    !   evaplg      Total particle loss rate due to evaporation for group
    !   evapls      Partial particle loss rate due to evaporation for element
    !   evappe      Particle production due to evaporation
    !   coreavg     Average total core mass in bin
    !   coresig     logarithm^2 of std dev of core distribution
    !   evdrop      Particle production of droplet number
    !   evcore      Particle production of core elements
    !   gasprod     Gas production term
    !   rlheat      Latent heating rate (per step) [deg_K/s]   
    !   ftoppart    Downward particle flux across top boundary of model
    !   fbotpart    Upward flux particle across bottom boundary of model
    !   pc_topbnd   Particle concentration assumed just above the top boundary
    !   pc_botbnd   Particle concentration assumed just below the bottom boundary
    !   cmf         Core mass fraction in a droplet 
    !   totevap     .true. if droplets are totally evaporating to CN
    !   too_small   .true. if cores are smaller than smallest CN
    !   too_big     .true. if cores are larger than largest CN
    !   nuc_small   .true. if cores are smaller than smallest nucleated CN
    !   rlprod      Latent heat production (per substep) (K/s)
    !
    real(kind=f), allocatable, dimension(:,:,:)   :: f_pcl        ! (NZ,NBIN,NELEM
    real(kind=f), allocatable, dimension(:,:)     :: f_gcl        ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)     :: f_d_gc       ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:)       :: f_d_t        ! (NZ)
    real(kind=f), allocatable, dimension(:,:)     :: f_dpc_sed    ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_pconmax    ! (NZ,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)   :: f_coaglg     ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)   :: f_coagpe     ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:,:)   :: f_rnuclg     ! (NBIN,NGROUP,NGROUP)
    real(kind=f), allocatable, dimension(:,:)     :: f_rnucpe     ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_rhompe     ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:,:)   :: f_pc_nucl    ! (NZ,NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_growpe     ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_evappe     ! (NBIN,NELEM)
    real(kind=f)                                  :: f_coreavg
    real(kind=f)                                  :: f_coresig
    real(kind=f)                                  :: f_evdrop
    real(kind=f), allocatable, dimension(:)       :: f_evcore     ! (NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_growlg     ! (NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)     :: f_evaplg     ! (NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:)       :: f_gasprod    ! (NGAS)
    real(kind=f), allocatable, dimension(:)       :: f_rlheat     ! (NZ)
    real(kind=f), allocatable, dimension(:,:)     :: f_ftoppart   ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_fbotpart   ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_pc_topbnd  ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_pc_botbnd  ! (NBIN,NELEM)
    real(kind=f), allocatable, dimension(:,:)     :: f_cmf        ! (NBIN,NGROUP)
    logical, allocatable, dimension(:,:)          :: f_totevap    ! (NBIN,NGROUP)
    logical                                       :: f_too_small
    logical                                       :: f_too_big
    logical                                       :: f_nuc_small
    real(kind=f)                                  :: f_rlprod

    !  Coagulation kernels and bin pair mapping
    !
    !   ckernel       Coagulation kernels [cm^3/s]          {setupckern}
    !
    real(kind=f), allocatable, dimension(:,:,:,:,:) :: f_ckernel ! (NZ,NBIN,NBIN,NGROUP,NGROUP)

    !  Particle fall velocities and diffusivities
    !
    !   bpm       Corrects for non-sphericity and non-continuum effects {setupvfall}
    !   vf        Fall velocities at layer endge                       {setupvfall}
    !   re        Reynolds' number based on <vfall>                     {setupvfall}
    !   dkz       Vert Brownian diffusion coef at layer boundary [z_units^2/s] {setupbdif}
    !   vd        Particle dry deposition velocity  [z_units/s]         {setupvdry}
    !
    real(kind=f), allocatable, dimension(:,:,:)     :: f_bpm        ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)     :: f_vf         ! (NZP1,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)     :: f_re         ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)     :: f_dkz        ! (NZP1,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)       :: f_vd         ! (NBIN,NGROUP)
    
    ! Atmospheric Structure
    !
    !  rhoa      Air density at layer mid-pt [g/x_units/y_units/z_units]  {initatm}
    !  rhoa_wet  Wet Air density averaged over grid box [g/x_units/y_units/z_units] {initatm}
    !  t         Air temperature at layer mid-pt [deg_K]                  {initatm}
    !  p         Atmospheric pressure at layer mid-pt [dyne/cm^2]         {initatm}
    !  pl        Atmospheric pressure at layer edge [dyne/cm^2]           {initatm}
    !  rmu       Air viscosity at layer mid-pt [g/cm/s]                   {initatm}
    !  thcond    Thermal conductivity of dry air [erg/cm/sec/deg_K]       {initatm}
    !  thcondnc  Adjusted thermal conductivity of dry air [erg/cm/sec/deg_K] {initatm}
    !  told      Temperature at beginning of time-step
    !  relhum    Hacked in relative humidity from hostmodel
    !  wtpct     Sulfate weight percent
    !
    real(kind=f), allocatable, dimension(:)     :: f_rhoa       ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_rhoa_wet   ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_t          ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_p          ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_pl         ! (NZP1)
    real(kind=f), allocatable, dimension(:)     :: f_rmu        ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_thcond     ! (NZ)
    real(kind=f), allocatable, dimension(:,:,:) :: f_thcondnc   ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:)     :: f_told       ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_relhum     ! (NZ)
    real(kind=f), allocatable, dimension(:)     :: f_wtpct      ! (NZ)

    ! Condensational growth parameters
    !
    ! NOTE: Some of these variables are used for storing intermediate values in
    ! the calculations. They may no longer be necessary, when the code is
    ! implemented as F90 and values as passed as parameters between subroutines.
    !
    !   diffus    Diffusivity of gas in air [cm^2/s]                  {setupgrow}
    !   rlhe      Latent heat of evaporation for gas [cm^2/s^2]       {setupgrow}
    !   rlhm      Latent heat of ice melting for gas [cm^2/s^2]       {setupgrow}
    !   pvapl     Saturation vapor pressure over water [dyne/cm^2]    {vaporp}   
    !   pvapi     Saturation vapor pressure over ice [dyne/cm^2]      {vaporp}   
    !   surfctwa  Surface tension of water-air interface              {setupgkern}
    !   surfctiw  Surface tension of water-ice interface              {setupgkern}
    !   surfctia  Surface tension of ice-air interface                {setupgkern}
    !   akelvin   Exponential arg. in curvature term for growth       {setupgkern}
    !   akelvini  Curvature term for ice                              {setupgkern}
    !   ft        Ventilation factor                                  {setupgkern}
    !   gro       Growth kernel [UNITS?]                              {setupgkern}
    !   gro1      Growth kernel conduction term [UNITS?]              {setupgkern}
    !   gro2      Growth kernel radiation term [UNITS?]               {setupgkern}
    !   supsatl   Supersaturation of vapor w.r.t. liquid water [dimless]
    !   supsati   Supersaturation of vapor w.r.t. ice [dimless]                  
    !   supsatlold Supersaturation (liquid) before time-step    {prestep}
    !   supsatiold Supersaturation (ice) before time-step    {prestep}
    !   scrit     Critical supersaturation for nucleation [dimless]   {setupnuc}
    !   radint    Incoming radiative intensity [erg/cm2/sr/s/um]
    !   partheat  Diffusional heating from particles (step) [K/s]
    !   dtpart    Delta particle temperature [K]
    !   phprod    Particle heating production (substep) [K/s]
    !
    real(kind=f), allocatable, dimension(:,:)    :: f_diffus     ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_rlhe       ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_rlhm       ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_pvapl      ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_pvapi      ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:)      :: f_surfctwa   ! (NZ)
    real(kind=f), allocatable, dimension(:)      :: f_surfctiw   ! (NZ)
    real(kind=f), allocatable, dimension(:)      :: f_surfctia   ! (NZ)
    real(kind=f), allocatable, dimension(:,:)    :: f_akelvin    ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_akelvini   ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:,:)  :: f_ft         ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)  :: f_gro        ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:,:)  :: f_gro1       ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)    :: f_gro2       ! (NZ,NGROUP)
    real(kind=f), allocatable, dimension(:,:)    :: f_supsatl    ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_supsati    ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_supsatlold ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:)    :: f_supsatiold ! (NZ,NGAS)
    real(kind=f), allocatable, dimension(:,:,:)  :: f_scrit      ! (NZ,NBIN,NGROUP)
    real(kind=f), allocatable, dimension(:,:)    :: f_radint     ! (NZ,NWAVE)
    real(kind=f), allocatable, dimension(:)      :: f_partheat   ! (NZ)
    real(kind=f), allocatable, dimension(:,:,:)  :: f_dtpart     ! (NZ,NBIN,NGROUP)
    real(kind=f)                                 :: f_phprod
   end type carmastate_type   
end module