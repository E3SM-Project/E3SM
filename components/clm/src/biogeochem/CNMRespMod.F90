module CNMRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding maintenance respiration routines for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varpar          , only : nlevgrnd
  use shr_const_mod       , only : SHR_CONST_TKFRZ
  use decompMod           , only : bounds_type
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use pftvarcon           , only : npcropmin
  use CNSharedParamsMod   , only : CNParamsShareInst
  use EcophysconType      , only : ecophyscon
  use SoilStateType       , only : soilstate_type
  use CanopyStateType     , only : canopystate_type
  use TemperatureType     , only : temperature_type
  use PhotosynthesisType  , only : photosyns_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PatchType           , only : pft                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNMResp
  public :: readCNMRespParams

  type, private :: CNMRespParamsType
     real(r8):: br        !base rate for maintenance respiration(gC/gN/s)
  end type CNMRespParamsType

  type(CNMRespParamsType),private ::  CNMRespParamsInst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNMRespParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read parameters
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNMRespParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='br_mr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNMRespParamsInst%br=tempr

  end subroutine readCNMRespParams

  !-----------------------------------------------------------------------
  ! FIX(SPM,032414) this shouldn't even be called with ED on.
  !
  subroutine CNMResp(bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       canopystate_vars, soilstate_vars, temperature_vars, photosyns_vars, &
       carbonflux_vars, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds          
    integer                  , intent(in)    :: num_soilc       ! number of soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_soilp       ! number of soil points in patch filter
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j ! indices
    integer :: fp    ! soil filter patch index
    integer :: fc    ! soil filter column index
    real(r8):: br    ! base rate (gC/gN/s)
    real(r8):: q10   ! temperature dependence
    real(r8):: tc    ! temperature correction, 2m air temp (unitless)
    real(r8):: tcsoi(bounds%begc:bounds%endc,nlevgrnd) ! temperature correction by soil layer (unitless)
    !-----------------------------------------------------------------------

    associate(                                                        &    
         ivt            =>    pft%itype                             , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         woody          =>    ecophyscon%woody                      , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)

         frac_veg_nosno =>    canopystate_vars%frac_veg_nosno_patch , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         laisun         =>    canopystate_vars%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index                  
         laisha         =>    canopystate_vars%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index                  

         rootfr         =>    soilstate_vars%rootfr_patch           , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer  (nlevgrnd)

         t_soisno       =>    temperature_vars%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         t_ref2m        =>    temperature_vars%t_ref2m_patch        , & ! Input:  [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)       

         lmrsun         =>    photosyns_vars%lmrsun_patch           , & ! Input:  [real(r8) (:)   ]  sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
         lmrsha         =>    photosyns_vars%lmrsha_patch           , & ! Input:  [real(r8) (:)   ]  shaded leaf maintenance respiration rate (umol CO2/m**2/s)

         leaf_mr        =>    carbonflux_vars%leaf_mr_patch         , & ! Output: [real(r8) (:)   ]                                                    
         froot_mr       =>    carbonflux_vars%froot_mr_patch        , & ! Output: [real(r8) (:)   ]                                                    
         livestem_mr    =>    carbonflux_vars%livestem_mr_patch     , & ! Output: [real(r8) (:)   ]                                                    
         livecroot_mr   =>    carbonflux_vars%livecroot_mr_patch    , & ! Output: [real(r8) (:)   ]                                                    
         grain_mr       =>    carbonflux_vars%grain_mr_patch        , & ! Output: [real(r8) (:)   ]                                                    

         frootn         =>    nitrogenstate_vars%frootn_patch       , & ! Input:  [real(r8) (:)   ]  (gN/m2) fine root N                               
         livestemn      =>    nitrogenstate_vars%livestemn_patch    , & ! Input:  [real(r8) (:)   ]  (gN/m2) live stem N                               
         livecrootn     =>    nitrogenstate_vars%livecrootn_patch   , & ! Input:  [real(r8) (:)   ]  (gN/m2) live coarse root N                        
         grainn         =>    nitrogenstate_vars%grainn_patch         & ! Output: [real(r8) (:)   ]  (kgN/m2) grain N
         )

      ! base rate for maintenance respiration is from:
      ! M. Ryan, 1991. Effects of climate change on plant respiration.
      ! Ecological Applications, 1(2), 157-167.
      ! Original expression is br = 0.0106 molC/(molN h)
      ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
      ! set constants
      br = CNMRespParamsInst%br

      ! Peter Thornton: 3/13/09 
      ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
      ! to improve seasonal cycle of atmospheric CO2 concentration in global
      ! simulatoins

      ! Set Q10 from CNSharedParamsMod
      Q10 = CNParamsShareInst%Q10

      ! column loop to calculate temperature factors in each soil layer
      do j=1,nlevgrnd
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! calculate temperature corrections for each soil layer, for use in
            ! estimating fine root maintenance respiration with depth
            tcsoi(c,j) = Q10**((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
         end do
      end do

      ! patch loop for leaves and live wood
      do fp = 1, num_soilp
         p = filter_soilp(fp)

         ! calculate maintenance respiration fluxes in
         ! gC/m2/s for each of the live plant tissues.
         ! Leaf and live wood MR

         tc = Q10**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)

         if (frac_veg_nosno(p) == 1) then

            leaf_mr(p) = lmrsun(p) * laisun(p) * 12.011e-6_r8 + &
                         lmrsha(p) * laisha(p) * 12.011e-6_r8

         else !nosno

            leaf_mr(p) = 0._r8

         end if

         if (woody(ivt(p)) == 1) then
            livestem_mr(p) = livestemn(p)*br*tc
            livecroot_mr(p) = livecrootn(p)*br*tc
         else if (ivt(p) >= npcropmin) then
            livestem_mr(p) = livestemn(p)*br*tc
            grain_mr(p) = grainn(p)*br*tc
         end if
      end do

      ! soil and patch loop for fine root

      do j = 1,nlevgrnd
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = pft%column(p)

            ! Fine root MR
            ! rootfr(j) sums to 1.0 over all soil layers, and
            ! describes the fraction of root mass that is in each
            ! layer.  This is used with the layer temperature correction
            ! to estimate the total fine root maintenance respiration as a
            ! function of temperature and N content.

            froot_mr(p) = froot_mr(p) + frootn(p)*br*tcsoi(c,j)*rootfr(p,j)
         end do
      end do

    end associate

  end subroutine CNMResp

end module CNMRespMod
