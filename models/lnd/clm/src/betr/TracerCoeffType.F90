 module TracerCoeffType
   use shr_kind_mod           , only : r8 => shr_kind_r8
   use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)  
   use decompMod              , only : bounds_type
   use PatchType              , only : pft                
   use ColumnType             , only : col                
   use LandunitType           , only : lun
   use landunit_varcon, only : istsoil, istcrop
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA:
  !
  !-------------------------------------------------------------------------------
  ! Column tracer phase conversion/transport parameters structure
  !-------------------------------------------------------------------------------
  type, public ::  TracerCoeff_type
    real(r8), pointer :: aqu2neutralcef_col        (:,:,:)      !aqueous tracer into neutral aqueous tracer
    real(r8), pointer :: aqu2bulkcef_mobile_col    (:,:,:)      !coefficient to convert bulk concentration into aqueous phase, (nlevsno+nlevtrc_soil)
    real(r8), pointer :: gas2bulkcef_mobile_col    (:,:,:)      !coefficient to convert bulk concentration into gaseous phase, (nlevsno+nlevlak+nlevtrc_soil)
    real(r8), pointer :: aqu2equilscef_col         (:,:,:)      !coefficient to convert solid phase (including ice) into aqueous phase
    real(r8), pointer :: henrycef_col              (:,:,:)      !henry's law constant
    real(r8), pointer :: bunsencef_col             (:,:,:)      !bunsen solubility
    real(r8), pointer :: tracer_diffusivity_air_col(:,:, :)     !diffusivity in the air
    real(r8), pointer :: aere_cond_col             (:,:)        !column level aerenchyma conductance (m/s)
    real(r8), pointer :: scal_aere_cond_col        (:,:)        !column level scaling factor for arenchyma or parenchyma transport
    real(r8), pointer :: diffgas_topsno_col        (:,:)        !gas diffusivity in top snow layer, this is not used currently
    real(r8), pointer :: diffgas_topsoi_col        (:,:)        !gas diffusivity in top soil layer, this is not used currently
    real(r8), pointer :: hmconductance_col         (:,:,:)      !geometrically weighted conductances (nlevsno+nlevtrc_soil)
    contains
      procedure, public  :: Init         
      procedure, private :: InitAllocate 
      procedure, private :: InitHistory  
      procedure, private :: InitCold      
  end type TracerCoeff_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, lbj, ubj, betrtracer_vars)
  
    use BeTRTracerType, only : BeTRTracer_Type  
    implicit none
    class(TracerCoeff_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars
    
    call this%InitAllocate(bounds, lbj, ubj, betrtracer_vars)
    call this%InitHistory(bounds, betrtracer_vars)
    call this%InitCold(bounds)

  end subroutine Init
  
  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, lbj, ubj, betrtracer_vars)

    use BeTRTracerType, only : BeTRTracer_Type  
    implicit none
    !
    ! !ARGUMENTS:
    class(TracerCoeff_type) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars    
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
 
    allocate(this%aqu2neutralcef_col         (begc:endc, lbj:ubj, 1:betrtracer_vars%ngwmobile_tracers));  this%aqu2neutralcef_col        (:,:,:) = nan
    allocate(this%aqu2bulkcef_mobile_col     (begc:endc, lbj:ubj, 1:betrtracer_vars%ngwmobile_tracers));  this%aqu2bulkcef_mobile_col    (:,:,:) = nan
    allocate(this%gas2bulkcef_mobile_col     (begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracers));  this%gas2bulkcef_mobile_col    (:,:,:) = nan
    allocate(this%henrycef_col               (begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracers));  this%henrycef_col              (:,:,:) = nan  
    allocate(this%bunsencef_col              (begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracers));  this%bunsencef_col             (:,:,:) = nan  
    allocate(this%tracer_diffusivity_air_col (begc:endc, lbj:ubj, 1:betrtracer_vars%nvolatile_tracers));  this%tracer_diffusivity_air_col(:,:,:) = nan  
    allocate(this%scal_aere_cond_col         (begc:endc,          1:betrtracer_vars%nvolatile_tracers));  this%scal_aere_cond_col        (:,:)   = nan
    allocate(this%aere_cond_col              (begc:endc,          1:betrtracer_vars%nvolatile_tracers));  this%aere_cond_col             (:,:)   = nan

    allocate(this%diffgas_topsno_col         (begc:endc,          1:betrtracer_vars%nvolatile_tracers));  this%diffgas_topsno_col        (:,:)   = nan
    allocate(this%diffgas_topsoi_col         (begc:endc,          1:betrtracer_vars%nvolatile_tracers));  this%diffgas_topsoi_col        (:,:)   = nan
    allocate(this%hmconductance_col          (begc:endc, lbj:ubj, 1:betrtracer_vars%ntracers));           this%hmconductance_col         (:,:,:) = nan
    allocate(this%aqu2equilscef_col          (begc:endc, lbj:ubj, 1:betrtracer_vars%nsolid_equil_tracers));this%aqu2equilscef_col        (:,:,:) = nan
  end subroutine InitAllocate
  
  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds, betrtracer_vars)
    !
    ! History fields initialization
    !
    ! !USES:
    !use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon     , only : spval
    use clm_varpar     , only : nlevsno
    use BeTRTracerType , only : BeTRTracer_Type
    use histFileMod    , only : hist_addfld1d, hist_addfld2d
    use histFileMod    , only : no_snow_normal, no_snow_zero
    
    !
    ! !ARGUMENTS:
    class(TracerCoeff_type)           :: this
    type(bounds_type)    , intent(in) :: bounds
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars        
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: jj, kk
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays
        
    !use the interface provided from CLM
    
    
    associate(                                                              &
      ntracers             => betrtracer_vars%ntracers                    , &
      ngwmobile_tracers    => betrtracer_vars%ngwmobile_tracers           , &
      nsolid_equil_tracers => betrtracer_vars%nsolid_equil_tracers        , &
      is_volatile          => betrtracer_vars%is_volatile                 , &
      volatileid           => betrtracer_vars%volatileid                  , &
      tracernames          => betrtracer_vars%tracernames                   &
    )
    begc=bounds%begc; endc=bounds%endc 
    do jj = 1, ntracers
      if(is_volatile(jj))then
        kk = volatileid(jj)
        this%scal_aere_cond_col(begc:endc, kk) = spval
        data1dptr => this%scal_aere_cond_col(begc:endc, kk)
        call hist_addfld1d (fname='SCAL_ARENCHYMA_'//tracernames(jj), units='none', &
         avgflag='A', long_name='scaling factor for tracer transport through arenchyma for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')

        this%aere_cond_col(begc:endc, kk) = spval
        data1dptr => this%aere_cond_col(begc:endc, kk)
        call hist_addfld1d (fname='ARENCHYMA_'//tracernames(jj), units='m/s', &
         avgflag='A', long_name='conductance for tracer transport through arenchyma for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')
         
        this%diffgas_topsoi_col(begc:endc, kk) = spval
        data1dptr => this%diffgas_topsoi_col(begc:endc, kk)
        call hist_addfld1d (fname='CDIFF_TOPSOI_'//tracernames(jj), units='none', &
         avgflag='A', long_name='gas diffusivity in top soil layer for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')

        this%gas2bulkcef_mobile_col(:,:,kk) = spval
        data2dptr => this%gas2bulkcef_mobile_col(:,:,kk)
        call hist_addfld2d (fname='CGAS2BULK_'//tracernames(jj), units='none', type2d='levtrc',  &
         avgflag='A', long_name='converting factor from gas to bulk phase for '//trim(tracernames(jj)), &
         ptr_col=data2dptr, default='inactive')        
      endif
      
      if(jj <= ngwmobile_tracers)then
        this%aqu2bulkcef_mobile_col(:,:,jj) = spval
        data2dptr => this%aqu2bulkcef_mobile_col(:,:,jj)
        call hist_addfld2d (fname='CAQU2BULK_vr_'//tracernames(jj), units='none', type2d='levtrc',  &
         avgflag='A', long_name='converting factor from aqeous to bulk phase for '//trim(tracernames(jj)), &
         ptr_col=data2dptr, default='inactive')
               
      endif
            
      this%hmconductance_col(:,:,jj) = spval
      data2dptr => this%hmconductance_col(:,:,jj)
      call hist_addfld2d (fname='HMCONDC_vr_'//tracernames(jj), units='none', type2d='levtrc',  &
       avgflag='A', long_name='bulk conductance for '//trim(tracernames(jj)), &
       ptr_col=data2dptr, default='inactive')              
    enddo     
        
    end associate    
  end subroutine InitHistory
  
  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use clm_varcon , only : spval
    !
    ! !ARGUMENTS:
    class(TracerCoeff_type) :: this
    type(bounds_type) , intent(in) :: bounds                                   
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l       ! index
    
    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (lun%ifspecial(l)) then
         this%aqu2neutralcef_col        (c,:,:) = spval
         this%aqu2bulkcef_mobile_col    (c,:,:) = spval
         this%gas2bulkcef_mobile_col    (c,:,:) = spval
         this%henrycef_col              (c,:,:) = spval
         this%bunsencef_col             (c,:,:) = spval
         this%tracer_diffusivity_air_col(c,:,:) = spval
         this%scal_aere_cond_col        (c,:)   = spval
         this%aere_cond_col             (c,:)   = spval
         this%diffgas_topsno_col        (c,:)   = spval
         this%diffgas_topsoi_col        (c,:)   = spval
         this%hmconductance_col         (c,:,:) = spval
       endif
       
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
         this%aqu2neutralcef_col        (c,:,:) = 0._r8
         this%aqu2bulkcef_mobile_col    (c,:,:) = 0._r8
         this%gas2bulkcef_mobile_col    (c,:,:) = 0._r8
         this%henrycef_col              (c,:,:) = 0._r8
         this%bunsencef_col             (c,:,:) = 0._r8
         this%tracer_diffusivity_air_col(c,:,:) = 0._r8
         this%scal_aere_cond_col        (c,:)   = 0._r8
         this%aere_cond_col             (c,:)   = 0._r8
         this%diffgas_topsno_col        (c,:)   = 0._r8
         this%diffgas_topsoi_col        (c,:)   = 0._r8
         this%hmconductance_col         (c,:,:) = 0._r8
       endif
    enddo
    
  end subroutine InitCold
  
  !------------------------------------------------------------------------

end module TracerCoeffType  
