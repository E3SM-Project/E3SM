module FrictionVelocityMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of the friction velocity, relation for potential
  ! temperature and humidity profiles of surface boundary layer.
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use clm_varcon   , only : spval
  use clm_varctl   , only: use_cn
  use LandunitType , only : lun                
  use ColumnType   , only : col
  use PatchType    , only : patch                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FrictionVelocity       ! Calculate friction velocity
  public :: MoninObukIni           ! Initialization of the Monin-Obukhov length
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: StabilityFunc1        ! Stability function for rib < 0.
  private :: StabilityFunc2        ! Stability function for rib < 0.

  type, public :: frictionvel_type

     ! Roughness length/resistance for friction velocity calculation

     real(r8), pointer, public :: forc_hgt_u_patch (:)   ! patch wind forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: forc_hgt_t_patch (:)   ! patch temperature forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: forc_hgt_q_patch (:)   ! patch specific humidity forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: u10_patch        (:)   ! patch 10-m wind (m/s) (for dust model)
     real(r8), pointer, public :: u10_clm_patch    (:)   ! patch 10-m wind (m/s) (for clm_map2gcell)
     real(r8), pointer, public :: va_patch         (:)   ! patch atmospheric wind speed plus convective velocity (m/s)
     real(r8), pointer, public :: vds_patch        (:)   ! patch deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
     real(r8), pointer, public :: fv_patch         (:)   ! patch friction velocity (m/s) (for dust model)
     real(r8), pointer, public :: rb1_patch        (:)   ! patch aerodynamical resistance (s/m) (for dry deposition of chemical tracers)
     real(r8), pointer, public :: ram1_patch       (:)   ! patch aerodynamical resistance (s/m)
     real(r8), pointer, public :: z0m_patch        (:)   ! patch momentum roughness length (m)
     real(r8), pointer, public :: z0mv_patch       (:)   ! patch roughness length over vegetation, momentum [m]
     real(r8), pointer, public :: z0hv_patch       (:)   ! patch roughness length over vegetation, sensible heat [m]
     real(r8), pointer, public :: z0qv_patch       (:)   ! patch roughness length over vegetation, latent heat [m]
     real(r8), pointer, public :: z0mg_col         (:)   ! col roughness length over ground, momentum  [m] 
     real(r8), pointer, public :: z0hg_col         (:)   ! col roughness length over ground, sensible heat [m]
     real(r8), pointer, public :: z0qg_col         (:)   ! col roughness length over ground, latent heat [m]

   contains

     ! Public procedures
     procedure, public  :: Init         
     procedure, public  :: Restart      

     ! Private procedures
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type frictionvel_type
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%forc_hgt_u_patch (begp:endp)) ; this%forc_hgt_u_patch (:)   = nan
    allocate(this%forc_hgt_t_patch (begp:endp)) ; this%forc_hgt_t_patch (:)   = nan
    allocate(this%forc_hgt_q_patch (begp:endp)) ; this%forc_hgt_q_patch (:)   = nan
    allocate(this%u10_patch        (begp:endp)) ; this%u10_patch        (:)   = nan
    allocate(this%u10_clm_patch    (begp:endp)) ; this%u10_clm_patch    (:)   = nan
    allocate(this%va_patch         (begp:endp)) ; this%va_patch         (:)   = nan
    allocate(this%vds_patch        (begp:endp)) ; this%vds_patch        (:)   = nan
    allocate(this%fv_patch         (begp:endp)) ; this%fv_patch         (:)   = nan
    allocate(this%rb1_patch        (begp:endp)) ; this%rb1_patch        (:)   = nan
    allocate(this%ram1_patch       (begp:endp)) ; this%ram1_patch       (:)   = nan
    allocate(this%z0m_patch        (begp:endp)) ; this%z0m_patch        (:)   = nan
    allocate(this%z0mv_patch       (begp:endp)) ; this%z0mv_patch       (:)   = nan
    allocate(this%z0hv_patch       (begp:endp)) ; this%z0hv_patch       (:)   = nan
    allocate(this%z0qv_patch       (begp:endp)) ; this%z0qv_patch       (:)   = nan
    allocate(this%z0mg_col         (begc:endc)) ; this%z0mg_col         (:)   = nan
    allocate(this%z0qg_col         (begc:endc)) ; this%z0qg_col         (:)   = nan
    allocate(this%z0hg_col         (begc:endc)) ; this%z0hg_col         (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%z0mg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0MG', units='m', &
         avgflag='A', long_name='roughness length over ground, momentum', &
         ptr_col=this%z0mg_col, default='inactive')

    this%z0hg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0HG', units='m', &
         avgflag='A', long_name='roughness length over ground, sensible heat', &
         ptr_col=this%z0hg_col, default='inactive')

    this%z0qg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0QG', units='m', &
         avgflag='A', long_name='roughness length over ground, latent heat', &
         ptr_col=this%z0qg_col, default='inactive')

    this%va_patch(begp:endp) = spval
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_patch=this%va_patch, default='inactive')

    this%u10_clm_patch(begp:endp) = spval
    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_patch=this%u10_clm_patch)

    if (use_cn) then
       this%u10_patch(begp:endp) = spval
       call hist_addfld1d (fname='U10_DUST', units='m/s', &
            avgflag='A', long_name='10-m wind for dust model', &
            ptr_patch=this%u10_patch, default='inactive')
    end if

    if (use_cn) then
       this%ram1_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAM1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%ram1_patch, default='inactive')
    end if

    if (use_cn) then
       this%fv_patch(begp:endp) = spval
       call hist_addfld1d (fname='FV', units='m/s', &
            avgflag='A', long_name='friction velocity for dust model', &
            ptr_patch=this%fv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0hv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0HV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, sensible heat', &
            ptr_patch=this%z0hv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0m_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0M', units='m', &
            avgflag='A', long_name='momentum roughness length', &
            ptr_patch=this%z0m_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0mv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0MV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, momentum', &
            ptr_patch=this%z0mv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0qv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0QV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, latent heat', &
            ptr_patch=this%z0qv_patch, default='inactive')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l                         ! indices
    !-----------------------------------------------------------------------

    ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
    ! since this is not initialized before first call to CNVegStructUpdate,
    ! and it is required to set the upper bound for canopy top height.
    ! Changed 3/21/08, KO: still needed but don't have sufficient information 
    ! to set this properly (e.g., patch-level displacement height and roughness 
    ! length). So leave at 30m.

    if (use_cn) then
       do p = bounds%begp, bounds%endp
          this%forc_hgt_u_patch(p) = 30._r8
       end do
    end if

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then !lake
          this%z0mg_col(c) = 0.0004_r8
       end if
    end do
   
  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='Z0MG', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground momentum roughness length', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%z0mg_col)

  end subroutine Restart

  !------------------------------------------------------------------------------
  subroutine FrictionVelocity(lbn, ubn, fn, filtern, &
       displa, z0m, z0h, z0q, &
       obu, iter, ur, um, ustar, &
       temp1, temp2, temp12m, temp22m, fm, frictionvel_inst, landunit_index)
    !
    ! !DESCRIPTION:
    ! Calculation of the friction velocity, relation for potential
    ! temperature and humidity profiles of surface boundary layer.
    ! The scheme is based on the work of Zeng et al. (1998):
    ! Intercomparison of bulk aerodynamic algorithms for the computation
    ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
    ! Vol. 11, 2628-2644.
    !
    ! !USES:
    use clm_varcon, only : vkc
    use clm_varctl, only : iulog
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: lbn, ubn                 ! pft/landunit array bounds
    integer  , intent(in)    :: fn                       ! number of filtered pft/landunit elements
    integer  , intent(in)    :: filtern(fn)              ! pft/landunit filter
    real(r8) , intent(in)    :: displa  ( lbn: )         ! displacement height (m) [lbn:ubn]
    real(r8) , intent(in)    :: z0m     ( lbn: )         ! roughness length over vegetation, momentum [m] [lbn:ubn]
    real(r8) , intent(in)    :: z0h     ( lbn: )         ! roughness length over vegetation, sensible heat [m] [lbn:ubn]
    real(r8) , intent(in)    :: z0q     ( lbn: )         ! roughness length over vegetation, latent heat [m] [lbn:ubn]
    real(r8) , intent(in)    :: obu     ( lbn: )         ! monin-obukhov length (m) [lbn:ubn]
    integer  , intent(in)    :: iter                     ! iteration number
    real(r8) , intent(in)    :: ur      ( lbn: )         ! wind speed at reference height [m/s] [lbn:ubn]
    real(r8) , intent(in)    :: um      ( lbn: )         ! wind speed including the stablity effect [m/s] [lbn:ubn]
    real(r8) , intent(out)   :: ustar   ( lbn: )         ! friction velocity [m/s] [lbn:ubn]
    real(r8) , intent(out)   :: temp1   ( lbn: )         ! relation for potential temperature profile [lbn:ubn]
    real(r8) , intent(out)   :: temp12m ( lbn: )         ! relation for potential temperature profile applied at 2-m [lbn:ubn]
    real(r8) , intent(out)   :: temp2   ( lbn: )         ! relation for specific humidity profile [lbn:ubn]
    real(r8) , intent(out)   :: temp22m ( lbn: )         ! relation for specific humidity profile applied at 2-m [lbn:ubn]
    real(r8) , intent(inout) :: fm      ( lbn: )         ! diagnose 10m wind (DUST only) [lbn:ubn]
    type(frictionvel_type) , intent(inout) :: frictionvel_inst
    logical  , intent(in), optional :: landunit_index   ! optional argument that defines landunit or pft level
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: zetam = 1.574_r8 ! transition point of flux-gradient relation (wind profile)
    real(r8), parameter :: zetat = 0.465_r8 ! transition point of flux-gradient relation (temp. profile)
    integer  :: f                           ! pft/landunit filter index
    integer  :: n                           ! pft/landunit index
    integer  :: g                           ! gridcell index
    integer  :: pp                          ! pfti,pftf index
    real(r8) :: zldis(lbn:ubn)              ! reference height "minus" zero displacement heght [m]
    real(r8) :: zeta(lbn:ubn)               ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: tmp1,tmp2,tmp3,tmp4         ! Used to diagnose the 10 meter wind
    real(r8) :: fmnew                       ! Used to diagnose the 10 meter wind
    real(r8) :: fm10                        ! Used to diagnose the 10 meter wind
    real(r8) :: zeta10                      ! Used to diagnose the 10 meter wind
    real(r8) :: vds_tmp                     ! Temporary for dry deposition velocity
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(displa)  == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(z0m)     == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(z0h)     == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(z0q)     == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(obu)     == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(ur)      == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(um)      == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(ustar)   == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(temp1)   == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(temp12m) == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(temp2)   == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(temp22m) == (/ubn/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fm)      == (/ubn/)), errMsg(__FILE__, __LINE__))

    associate(                                                   & 
         pfti             => lun%patchi                          , & ! Input:  [integer  (:) ] beginning pfti index for landunit         
         pftf             => lun%patchf                          , & ! Input:  [integer  (:) ] final pft index for landunit              
         
         forc_hgt_u_patch => frictionvel_inst%forc_hgt_u_patch , & ! Input:  [real(r8) (:) ] observational height of wind at pft level [m]
         forc_hgt_t_patch => frictionvel_inst%forc_hgt_t_patch , & ! Input:  [real(r8) (:) ] observational height of temperature at pft level [m]
         forc_hgt_q_patch => frictionvel_inst%forc_hgt_q_patch , & ! Input:  [real(r8) (:) ] observational height of specific humidity at pft level [m]
         vds              => frictionvel_inst%vds_patch        , & ! Output: [real(r8) (:) ] dry deposition velocity term (m/s) (for SO4 NH4NO3)
         u10              => frictionvel_inst%u10_patch        , & ! Output: [real(r8) (:) ] 10-m wind (m/s) (for dust model)        
         u10_clm          => frictionvel_inst%u10_clm_patch    , & ! Output: [real(r8) (:) ] 10-m wind (m/s)                         
         va               => frictionvel_inst%va_patch         , & ! Output: [real(r8) (:) ] atmospheric wind speed plus convective velocity (m/s)
         fv               => frictionvel_inst%fv_patch           & ! Output: [real(r8) (:) ] friction velocity (m/s) (for dust model)
         )

      ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

      do f = 1, fn
         n = filtern(f)
         if (present(landunit_index)) then
            g = lun%gridcell(n) 
         else
            g = patch%gridcell(n)  
         end if

         ! Wind profile

         if (present(landunit_index)) then
            zldis(n) = forc_hgt_u_patch(pfti(n))-displa(n)
         else
            zldis(n) = forc_hgt_u_patch(n)-displa(n)
         end if
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetam) then
            ustar(n) = vkc*um(n)/(log(-zetam*obu(n)/z0m(n))&
                 - StabilityFunc1(-zetam) &
                 + StabilityFunc1(z0m(n)/obu(n)) &
                 + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8))
         else if (zeta(n) < 0._r8) then
            ustar(n) = vkc*um(n)/(log(zldis(n)/z0m(n))&
                 - StabilityFunc1(zeta(n))&
                 + StabilityFunc1(z0m(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            ustar(n) = vkc*um(n)/(log(zldis(n)/z0m(n)) + 5._r8*zeta(n) -5._r8*z0m(n)/obu(n))
         else
            ustar(n) = vkc*um(n)/(log(obu(n)/z0m(n))+5._r8-5._r8*z0m(n)/obu(n) &
                 +(5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if

         if (zeta(n) < 0._r8) then
            vds_tmp = 2.e-3_r8*ustar(n) * ( 1._r8 + (300._r8/(-obu(n)))**0.666_r8)
         else
            vds_tmp = 2.e-3_r8*ustar(n)
         endif

         if (present(landunit_index)) then
            do pp = pfti(n),pftf(n)
               vds(pp) = vds_tmp
            end do
         else
            vds(n) = vds_tmp
         end if

         ! Calculate a 10-m wind (10m + z0m + d)
         ! For now, this will not be the same as the 10-m wind calculated for the dust 
         ! model because the CLM stability functions are used here, not the LSM stability
         ! functions used in the dust model. We will eventually change the dust model to be 
         ! consistent with the following formulation.
         ! Note that the 10-m wind calculated this way could actually be larger than the
         ! atmospheric forcing wind because 1) this includes the convective velocity, 2)
         ! this includes the 1 m/s minimum wind threshold

         ! If forcing height is less than or equal to 10m, then set 10-m wind to um
         if (present(landunit_index)) then
            do pp = pfti(n),pftf(n)
               if (zldis(n)-z0m(n) <= 10._r8) then
                  u10_clm(pp) = um(n)
               else
                  if (zeta(n) < -zetam) then
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))      &
                          - StabilityFunc1(-zetam)                              &
                          + StabilityFunc1((10._r8+z0m(n))/obu(n))              &
                          + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
                  else if (zeta(n) < 0._r8) then
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                          - StabilityFunc1(zeta(n))                             &
                          + StabilityFunc1((10._r8+z0m(n))/obu(n))) )
                  else if (zeta(n) <=  1._r8) then
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                          + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
                  else
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))             &
                          + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                &
                          + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )

                  end if
               end if
               va(pp) = um(n)
            end do
         else
            if (zldis(n)-z0m(n) <= 10._r8) then
               u10_clm(n) = um(n)
            else
               if (zeta(n) < -zetam) then
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))         &
                       - StabilityFunc1(-zetam)                                 &
                       + StabilityFunc1((10._r8+z0m(n))/obu(n))                 &
                       + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
               else if (zeta(n) < 0._r8) then
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                       - StabilityFunc1(zeta(n))                                &
                       + StabilityFunc1((10._r8+z0m(n))/obu(n))) )
               else if (zeta(n) <=  1._r8) then
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                       + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
               else
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))                &
                       + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                   &
                       + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )
               end if
            end if
            va(n) = um(n)
         end if

         ! Temperature profile

         if (present(landunit_index)) then
            zldis(n) = forc_hgt_t_patch(pfti(n))-displa(n)
         else
            zldis(n) = forc_hgt_t_patch(n)-displa(n)
         end if
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp1(n) = vkc/(log(-zetat*obu(n)/z0h(n))&
                 - StabilityFunc2(-zetat) &
                 + StabilityFunc2(z0h(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp1(n) = vkc/(log(zldis(n)/z0h(n)) &
                 - StabilityFunc2(zeta(n)) &
                 + StabilityFunc2(z0h(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp1(n) = vkc/(log(zldis(n)/z0h(n)) + 5._r8*zeta(n) - 5._r8*z0h(n)/obu(n))
         else
            temp1(n) = vkc/(log(obu(n)/z0h(n)) + 5._r8 - 5._r8*z0h(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if

         ! Humidity profile

         if (present(landunit_index)) then
            if (forc_hgt_q_patch(pfti(n)) == forc_hgt_t_patch(pfti(n)) .and. z0q(n) == z0h(n)) then
               temp2(n) = temp1(n)
            else
               zldis(n) = forc_hgt_q_patch(pfti(n))-displa(n)
               zeta(n) = zldis(n)/obu(n)
               if (zeta(n) < -zetat) then
                  temp2(n) = vkc/(log(-zetat*obu(n)/z0q(n)) &
                       - StabilityFunc2(-zetat) &
                       + StabilityFunc2(z0q(n)/obu(n)) &
                       + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
               else if (zeta(n) < 0._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) &
                       - StabilityFunc2(zeta(n)) &
                       + StabilityFunc2(z0q(n)/obu(n)))
               else if (zeta(n) <=  1._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
               else
                  temp2(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                       + (5._r8*log(zeta(n))+zeta(n)-1._r8))
               end if
            end if
         else
            if (forc_hgt_q_patch(n) == forc_hgt_t_patch(n) .and. z0q(n) == z0h(n)) then
               temp2(n) = temp1(n)
            else
               zldis(n) = forc_hgt_q_patch(n)-displa(n)
               zeta(n) = zldis(n)/obu(n)
               if (zeta(n) < -zetat) then
                  temp2(n) = vkc/(log(-zetat*obu(n)/z0q(n)) &
                       - StabilityFunc2(-zetat) &
                       + StabilityFunc2(z0q(n)/obu(n)) &
                       + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
               else if (zeta(n) < 0._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) &
                       - StabilityFunc2(zeta(n)) &
                       + StabilityFunc2(z0q(n)/obu(n)))
               else if (zeta(n) <=  1._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
               else
                  temp2(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                       + (5._r8*log(zeta(n))+zeta(n)-1._r8))
               end if
            endif
         endif

         ! Temperature profile applied at 2-m

         zldis(n) = 2.0_r8 + z0h(n)
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp12m(n) = vkc/(log(-zetat*obu(n)/z0h(n))&
                 - StabilityFunc2(-zetat) &
                 + StabilityFunc2(z0h(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp12m(n) = vkc/(log(zldis(n)/z0h(n)) &
                 - StabilityFunc2(zeta(n))  &
                 + StabilityFunc2(z0h(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp12m(n) = vkc/(log(zldis(n)/z0h(n)) + 5._r8*zeta(n) - 5._r8*z0h(n)/obu(n))
         else
            temp12m(n) = vkc/(log(obu(n)/z0h(n)) + 5._r8 - 5._r8*z0h(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if

         ! Humidity profile applied at 2-m

         if (z0q(n) == z0h(n)) then
            temp22m(n) = temp12m(n)
         else
            zldis(n) = 2.0_r8 + z0q(n)
            zeta(n) = zldis(n)/obu(n)
            if (zeta(n) < -zetat) then
               temp22m(n) = vkc/(log(-zetat*obu(n)/z0q(n)) - &
                    StabilityFunc2(-zetat) + StabilityFunc2(z0q(n)/obu(n)) &
                    + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
            else if (zeta(n) < 0._r8) then
               temp22m(n) = vkc/(log(zldis(n)/z0q(n)) - &
                    StabilityFunc2(zeta(n))+StabilityFunc2(z0q(n)/obu(n)))
            else if (zeta(n) <=  1._r8) then
               temp22m(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
            else
               temp22m(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                    + (5._r8*log(zeta(n))+zeta(n)-1._r8))
            end if
         end if

         ! diagnose 10-m wind for dust model (dstmbl.F)
         ! Notes from C. Zender's dst.F:
         ! According to Bon96 p. 62, the displacement height d (here displa) is
         ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
         ! Therefore d <= 0.034*z1 and may safely be neglected.
         ! Code from LSM routine SurfaceTemperature was used to obtain u10

         if (present(landunit_index)) then
            zldis(n) = forc_hgt_u_patch(pfti(n))-displa(n)
         else
            zldis(n) = forc_hgt_u_patch(n)-displa(n)
         end if
         zeta(n) = zldis(n)/obu(n)
         if (min(zeta(n), 1._r8) < 0._r8) then
            tmp1 = (1._r8 - 16._r8*min(zeta(n),1._r8))**0.25_r8
            tmp2 = log((1._r8+tmp1*tmp1)/2._r8)
            tmp3 = log((1._r8+tmp1)/2._r8)
            fmnew = 2._r8*tmp3 + tmp2 - 2._r8*atan(tmp1) + 1.5707963_r8
         else
            fmnew = -5._r8*min(zeta(n),1._r8)
         endif
         if (iter == 1) then
            fm(n) = fmnew
         else
            fm(n) = 0.5_r8 * (fm(n)+fmnew)
         end if
         zeta10 = min(10._r8/obu(n), 1._r8)
         if (zeta(n) == 0._r8) zeta10 = 0._r8
         if (zeta10 < 0._r8) then
            tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
            tmp2 = log((1.0_r8 + tmp1*tmp1)/2.0_r8)
            tmp3 = log((1.0_r8 + tmp1)/2.0_r8)
            fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*atan(tmp1) + 1.5707963_r8
         else                ! not stable
            fm10 = -5.0_r8 * zeta10
         end if
         if (present(landunit_index)) then
            tmp4 = log( max( 1.0_r8, forc_hgt_u_patch(pfti(n)) / 10._r8) )
         else 
            tmp4 = log( max( 1.0_r8, forc_hgt_u_patch(n) / 10._r8) )
         end if
         if (present(landunit_index)) then
            do pp = pfti(n),pftf(n)
               u10(pp) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
               fv(pp)  = ustar(n)
            end do
         else
            u10(n) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
            fv(n)  = ustar(n)
         end if

      end do

    end associate
  end subroutine FrictionVelocity

  !------------------------------------------------------------------------------
  real(r8) function StabilityFunc1(zeta)
    !
    ! !DESCRIPTION:
    ! Stability function for rib < 0.
    !
    ! !USES:
    use shr_const_mod, only: SHR_CONST_PI
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    !
    ! !LOCAL VARIABLES:
    real(r8) :: chik, chik2
    !------------------------------------------------------------------------------

    chik2 = sqrt(1._r8-16._r8*zeta)
    chik = sqrt(chik2)
    StabilityFunc1 = 2._r8*log((1._r8+chik)*0.5_r8) &
         + log((1._r8+chik2)*0.5_r8)-2._r8*atan(chik)+SHR_CONST_PI*0.5_r8

  end function StabilityFunc1

  !------------------------------------------------------------------------------
  real(r8) function StabilityFunc2(zeta)
    !
    ! !DESCRIPTION:
    ! Stability function for rib < 0.
    !
    ! !USES:
    use shr_const_mod, only: SHR_CONST_PI
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    !
    ! !LOCAL VARIABLES:
    real(r8) :: chik2
    !------------------------------------------------------------------------------

    chik2 = sqrt(1._r8-16._r8*zeta)
    StabilityFunc2 = 2._r8*log((1._r8+chik2)*0.5_r8)

  end function StabilityFunc2

  !-----------------------------------------------------------------------
  subroutine MoninObukIni (ur, thv, dthv, zldis, z0m, um, obu)
    !
    ! !DESCRIPTION:
    ! Initialization of the Monin-Obukhov length.
    ! The scheme is based on the work of Zeng et al. (1998):
    ! Intercomparison of bulk aerodynamic algorithms for the computation
    ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
    ! Vol. 11, 2628-2644.
    !
    ! !USES:
    use clm_varcon, only : grav
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in)  :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in)  :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in)  :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in)  :: z0m   ! roughness length, momentum [m]
    real(r8), intent(out) :: um    ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu   ! monin-obukhov length (m)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wc    ! convective velocity [m/s]
    real(r8) :: rib   ! bulk Richardson number
    real(r8) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: ustar ! friction velocity [m/s]
    !-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06_r8
    wc=0.5_r8
    if (dthv >= 0._r8) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
       zeta = min(2._r8,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

  end subroutine MoninObukIni

end module FrictionVelocityMod
