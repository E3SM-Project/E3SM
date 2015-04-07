module SoilBiogeochemStateType

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevsoifl, nlevsoi, crop_prog
  use clm_varpar     , only : ndecomp_cascade_transitions, nlevdecomp, nlevdecomp_full, more_vertlayers  
  use clm_varcon     , only : spval, ispval, c14ratio, grlnd
  use landunit_varcon, only : istsoil, istcrop
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, crop_prog 
  use clm_varctl     , only : use_vertsoilc, use_cn 
  use clm_varctl     , only : iulog
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilbiogeochem_state_type

     real(r8) , pointer :: leaf_prof_patch             (:,:)   ! (1/m) profile of leaves (vertical profiles for calculating fluxes)
     real(r8) , pointer :: froot_prof_patch            (:,:)   ! (1/m) profile of fine roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: croot_prof_patch            (:,:)   ! (1/m) profile of coarse roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: stem_prof_patch             (:,:)   ! (1/m) profile of stems (vertical profiles for calculating fluxes)
     real(r8) , pointer :: fpi_vr_col                  (:,:)   ! (no units) fraction of potential immobilization 
     real(r8) , pointer :: fpi_col                     (:)     ! (no units) fraction of potential immobilization 
     real(r8),  pointer :: fpg_col                     (:)     ! (no units) fraction of potential gpp 
     real(r8) , pointer :: rf_decomp_cascade_col       (:,:,:) ! (frac) respired fraction in decomposition step 
     real(r8) , pointer :: pathfrac_decomp_cascade_col (:,:,:) ! (frac) what fraction of C leaving a given pool passes through a given transition 
     real(r8) , pointer :: nfixation_prof_col          (:,:)   ! (1/m) profile for N fixation additions 
     real(r8) , pointer :: ndep_prof_col               (:,:)   ! (1/m) profile for N fixation additions 
     real(r8) , pointer :: som_adv_coef_col            (:,:)   ! (m2/s) SOM advective flux 
     real(r8) , pointer :: som_diffus_coef_col         (:,:)   ! (m2/s) SOM diffusivity due to bio/cryo-turbation 
     real(r8) , pointer :: plant_ndemand_col           (:)     ! column-level plant N demand

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type soilbiogeochem_state_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    if (use_cn) then
       call this%InitHistory ( bounds )
    end if
    call this%InitCold ( bounds ) 

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
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%leaf_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%leaf_prof_patch     (:,:) = spval
    allocate(this%froot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%froot_prof_patch    (:,:) = spval
    allocate(this%croot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%croot_prof_patch    (:,:) = spval
    allocate(this%stem_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%stem_prof_patch     (:,:) = spval
    allocate(this%fpi_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%fpi_vr_col          (:,:) = nan
    allocate(this%fpi_col             (begc:endc))                   ; this%fpi_col             (:)   = nan
    allocate(this%fpg_col             (begc:endc))                   ; this%fpg_col             (:)   = nan
    allocate(this%nfixation_prof_col  (begc:endc,1:nlevdecomp_full)) ; this%nfixation_prof_col  (:,:) = spval
    allocate(this%ndep_prof_col       (begc:endc,1:nlevdecomp_full)) ; this%ndep_prof_col       (:,:) = spval
    allocate(this%som_adv_coef_col    (begc:endc,1:nlevdecomp_full)) ; this%som_adv_coef_col    (:,:) = spval
    allocate(this%som_diffus_coef_col (begc:endc,1:nlevdecomp_full)) ; this%som_diffus_coef_col (:,:) = spval
    allocate(this%plant_ndemand_col   (begc:endc))                   ; this%plant_ndemand_col   (:)   = nan

    allocate(this%rf_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); 
    this%rf_decomp_cascade_col(:,:,:) = nan

    allocate(this%pathfrac_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));     
    this%pathfrac_decomp_cascade_col(:,:,:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp, no_snow_normal
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    character(8)      :: vr_suffix
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%croot_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='CROOT_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from coarse roots', &
         ptr_patch=this%croot_prof_patch, default='inactive')

    this%froot_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='FROOT_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from fine roots', &
         ptr_patch=this%froot_prof_patch, default='inactive')

    this%leaf_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='LEAF_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from leaves', &
         ptr_patch=this%leaf_prof_patch, default='inactive')

    this%stem_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='STEM_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from stems', &
         ptr_patch=this%stem_prof_patch, default='inactive')

    this%nfixation_prof_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='NFIXATION_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for biological N fixation', &
         ptr_col=this%nfixation_prof_col, default='inactive')

    this%ndep_prof_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='NDEP_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for atmospheric N  deposition', &
         ptr_col=this%ndep_prof_col, default='inactive')

    this%som_adv_coef_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SOM_ADV_COEF', units='m/s',  type2d='levdcmp', &
         avgflag='A', long_name='advection term for vertical SOM translocation', &
         ptr_col=this%som_adv_coef_col, default='inactive')

    this%som_diffus_coef_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SOM_DIFFUS_COEF', units='m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='diffusion coefficient for vertical SOM translocation', &
         ptr_col=this%som_diffus_coef_col, default='inactive')

    if ( nlevdecomp_full > 1 ) then
       this%fpi_col(begc:endc) = spval
       call hist_addfld1d (fname='FPI', units='proportion', &
            avgflag='A', long_name='fraction of potential immobilization', &
            ptr_col=this%fpi_col)
    endif

    this%fpg_col(begc:endc) = spval
    call hist_addfld1d (fname='FPG', units='proportion', &
         avgflag='A', long_name='fraction of potential gpp', &
         ptr_col=this%fpg_col)

    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif
    this%fpi_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='FPI'//trim(vr_suffix), units='proportion', type2d='levdcmp', & 
         avgflag='A', long_name='fraction of potential immobilization', &
         ptr_col=this%fpi_vr_col)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use clm_varctl , only : fsurdat
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds   
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    type(file_desc_t)     :: ncid                     ! netcdf id
    logical               :: readvar 
    character(len=256)    :: locfn                    ! local filename
    integer               :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    ! --------------------------------------------------------------------
    ! Open surface dataset
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_inqdlen(ncid,dimid,nlevsoifl,name='nlevsoi')
    if ( .not. more_vertlayers )then
       if ( nlevsoifl /= nlevsoi )then
          call endrun(msg=' ERROR: Number of soil layers on file does NOT match the number being used'//&
               errMsg(__FILE__, __LINE__))
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if

    ! --------------------------------------------------------------------
    ! Initialize terms needed for dust model
    ! --------------------------------------------------------------------
       
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          this%fpi_col (c) = spval
          this%fpg_col (c) = spval
          do j = 1,nlevdecomp_full
             this%fpi_vr_col(c,j) = spval
          end do
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          ! initialize fpi_vr so that levels below nlevsoi are not nans
          this%fpi_vr_col(c,1:nlevdecomp_full)          = 0._r8 
          this%som_adv_coef_col(c,1:nlevdecomp_full)    = 0._r8 
          this%som_diffus_coef_col(c,1:nlevdecomp_full) = 0._r8 

          ! initialize the profiles for converting to vertically resolved carbon pools
          this%nfixation_prof_col(c,1:nlevdecomp_full)  = 0._r8 
          this%ndep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
       end if
    end do

  end subroutine initCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer, pointer  :: temp1d(:)  ! temporary
    integer           :: p,j,c,i    ! indices
    logical           :: readvar    ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    !-----------------------------------------------------------------------
  
    if (use_vertsoilc) then
       ptr2d => this%fpi_vr_col
       call restartvar(ncid=ncid, flag=flag, varname='fpi_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='fraction of potential immobilization',  units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%fpi_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
       call restartvar(ncid=ncid, flag=flag, varname='fpi', xtype=ncd_double,  &
            dim1name='column', &
            long_name='fraction of potential immobilization',  units='unitless', &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    if (use_vertsoilc) then
       ptr2d => this%som_adv_coef_col
       call restartvar(ncid=ncid, flag=flag, varname='som_adv_coef_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='SOM advective flux', units='m/s', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if
    
    if (use_vertsoilc) then
       ptr2d => this%som_diffus_coef_col
       call restartvar(ncid=ncid, flag=flag, varname='som_diffus_coef_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='SOM diffusivity due to bio/cryo-turbation',  units='m^2/s', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fpg', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpg_col) 

  end subroutine Restart

end module SoilBiogeochemStateType
