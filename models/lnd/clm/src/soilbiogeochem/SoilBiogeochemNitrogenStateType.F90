module SoilBiogeochemNitrogenStateType

#include "shr_assert.h"

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use decompMod                          , only : bounds_type
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc 
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon                         , only : spval, dzsoi_decomp, zisoi
  use clm_varctl                         , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varctl                         , only : iulog, override_bgc_restart_mismatch_dump, spinup_state
  use landunit_varcon                    , only : istcrop, istsoil 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: soilbiogeochem_nitrogenstate_type

     real(r8), pointer :: decomp_npools_vr_col         (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: sminn_vr_col                 (:,:)   ! col (gN/m3) vertically-resolved soil mineral N
     real(r8), pointer :: ntrunc_vr_col                (:,:)   ! col (gN/m3) vertically-resolved column-level sink for N truncation

     ! nitrif_denitrif
     real(r8), pointer :: smin_no3_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_no3_col                 (:)     ! col (gN/m2) soil mineral NO3 pool
     real(r8), pointer :: smin_nh4_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4_col                 (:)     ! col (gN/m2) soil mineral NH4 pool

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: decomp_npools_col            (:,:)   ! col (gN/m2)  decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: decomp_npools_1m_col         (:,:)   ! col (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
     real(r8), pointer :: sminn_col                    (:)     ! col (gN/m2) soil mineral N
     real(r8), pointer :: ntrunc_col                   (:)     ! col (gN/m2) column-level sink for N truncation
     real(r8), pointer :: cwdn_col                     (:)     ! col (gN/m2) Diagnostic: coarse woody debris N
     real(r8), pointer :: totlitn_col                  (:)     ! col (gN/m2) total litter nitrogen
     real(r8), pointer :: totsomn_col                  (:)     ! col (gN/m2) total soil organic matter nitrogen
     real(r8), pointer :: totlitn_1m_col               (:)     ! col (gN/m2) total litter nitrogen to 1 meter
     real(r8), pointer :: totsomn_1m_col               (:)     ! col (gN/m2) total soil organic matter nitrogen to 1 meter

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: Summary
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type soilbiogeochem_nitrogenstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)

    class(soilbiogeochem_nitrogenstate_type)         :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: decomp_cpools_vr_col (bounds%begc:, 1:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_col    (bounds%begc:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_1m_col (bounds%begc:, 1:)

    call this%InitAllocate (bounds )

    call this%InitHistory (bounds)

    call this%InitCold ( bounds, & 
         decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begc,endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    allocate(this%sminn_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%sminn_vr_col         (:,:) = nan
    allocate(this%ntrunc_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%ntrunc_vr_col        (:,:) = nan
    allocate(this%smin_no3_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_vr_col      (:,:) = nan
    allocate(this%smin_nh4_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_vr_col      (:,:) = nan
    allocate(this%smin_no3_col         (begc:endc))                   ; this%smin_no3_col         (:)   = nan
    allocate(this%smin_nh4_col         (begc:endc))                   ; this%smin_nh4_col         (:)   = nan
    allocate(this%cwdn_col             (begc:endc))                   ; this%cwdn_col             (:)   = nan
    allocate(this%sminn_col            (begc:endc))                   ; this%sminn_col            (:)   = nan
    allocate(this%ntrunc_col           (begc:endc))                   ; this%ntrunc_col           (:)   = nan
    allocate(this%totlitn_col          (begc:endc))                   ; this%totlitn_col          (:)   = nan
    allocate(this%totsomn_col          (begc:endc))                   ; this%totsomn_col          (:)   = nan
    allocate(this%totlitn_1m_col       (begc:endc))                   ; this%totlitn_1m_col       (:)   = nan
    allocate(this%totsomn_1m_col       (begc:endc))                   ; this%totsomn_1m_col       (:)   = nan
    allocate(this%decomp_npools_col    (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_col    (:,:) = nan
    allocate(this%decomp_npools_1m_col (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_1m_col (:,:) = nan

    allocate(this%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%decomp_npools_vr_col(:,:,:)= nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full,crop_prog, nlevgrnd
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    character(8)      :: vr_suffix
    integer           :: begc,endc
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    if ( nlevdecomp_full > 1 ) then
       this%decomp_npools_vr_col(begc:endc,:,:) = spval
       this%decomp_npools_1m_col(begc:endc,:) = spval
    end if
    this%decomp_npools_col(begc:endc,:) = spval
    do l  = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => this%decomp_npools_vr_col(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif

       data1dptr => this%decomp_npools_col(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N'
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N'
       call hist_addfld1d (fname=fieldname, units='gN/m^2', &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

       if ( nlevdecomp_full > 1 ) then
          data1dptr => this%decomp_npools_1m_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N to 1 meter'
          call hist_addfld1d (fname=fieldname, units='gN/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif
    end do


    if ( nlevdecomp_full > 1 ) then

       this%sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINN', units='gN/m^2', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_col)

       this%totlitn_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTLITN_1m', units='gN/m^2', &
            avgflag='A', long_name='total litter N to 1 meter', &
            ptr_col=this%totlitn_1m_col)

       this%totsomn_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMN_1m', units='gN/m^2', &
            avgflag='A', long_name='total soil organic matter N to 1 meter', &
            ptr_col=this%totsomn_1m_col)
    endif

    this%ntrunc_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_NTRUNC', units='gN/m^2',  &
         avgflag='A', long_name='column-level sink for N truncation', &
         ptr_col=this%ntrunc_col)

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    if (use_nitrif_denitrif) then
       this%smin_no3_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NO3 (vert. res.)', &
            ptr_col=this%smin_no3_vr_col)

       this%smin_nh4_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NH4'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NH4 (vert. res.)', &
            ptr_col=this%smin_nh4_vr_col)

       if ( nlevdecomp_full > 1 ) then
          this%smin_no3_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NO3', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NO3', &
               ptr_col=this%smin_no3_col)

          this%smin_nh4_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NH4', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NH4', &
               ptr_col=this%smin_nh4_col)
       endif

       this%sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_vr_col, default = 'inactive')
    else
       this%sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_vr_col)
    end if

    this%totlitn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTLITN', units='gN/m^2', &
         avgflag='A', long_name='total litter N', &
         ptr_col=this%totlitn_col)

    this%totsomn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTSOMN', units='gN/m^2', &
         avgflag='A', long_name='total soil organic matter N', &
         ptr_col=this%totsomn_col)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use clm_varpar     , only : crop_prog
    use decompMod      , only : bounds_type
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenstate_type)      :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: decomp_cpools_vr_col(bounds%begc:,:,:)
    real(r8)          , intent(in) :: decomp_cpools_col(bounds%begc:,:)
    real(r8)          , intent(in) :: decomp_cpools_1m_col(bounds%begc:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,g,l,c,j,k                              ! indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col   (bounds%endc-bounds%begc+1) ! special landunit filter - columns
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(decomp_cpools_col)    == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_1m_col) == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)), errMsg(__FILE__, __LINE__))

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          ! column nitrogen state variables
          this%ntrunc_col(c) = 0._r8
          this%sminn_col(c) = 0._r8
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                this%decomp_npools_vr_col(c,j,k) = decomp_cpools_vr_col(c,j,k) / decomp_cascade_con%initial_cn_ratio(k)
             end do
             this%sminn_vr_col(c,j) = 0._r8
             this%ntrunc_vr_col(c,j) = 0._r8
          end do
          if ( nlevdecomp > 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   this%decomp_npools_vr_col(c,j,k) = 0._r8
                end do
                this%sminn_vr_col(c,j) = 0._r8
                this%ntrunc_vr_col(c,j) = 0._r8
             end do
          end if
          do k = 1, ndecomp_pools
             this%decomp_npools_col(c,k)    = decomp_cpools_col(c,k)    / decomp_cascade_con%initial_cn_ratio(k)
             this%decomp_npools_1m_col(c,k) = decomp_cpools_1m_col(c,k) / decomp_cascade_con%initial_cn_ratio(k)
          end do

          if (use_nitrif_denitrif) then
             do j = 1, nlevdecomp_full
                this%smin_nh4_vr_col(c,j) = 0._r8
                this%smin_no3_vr_col(c,j) = 0._r8
             end do
             this%smin_nh4_col(c) = 0._r8
             this%smin_no3_col(c) = 0._r8
          end if
          this%totlitn_col(c)    = 0._r8
          this%totsomn_col(c)    = 0._r8
          this%totlitn_1m_col(c) = 0._r8
          this%totsomn_1m_col(c) = 0._r8
          this%cwdn_col(c)       = 0._r8

       end if
    end do

    ! initialize fields for special filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    call this%SetValues (num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod      , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_time_manager    , only : is_restart, get_nstep
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m          ! multiplier for the exit_spinup code
    real(r8), pointer  :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    integer            :: itemp      ! temporary 
    integer , pointer  :: iptemp(:)  ! pointer to memory to be allocated
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    ! sminn
    if (use_vertsoilc) then
       ptr2d => this%sminn_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="sminn_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%sminn_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="sminn", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR::'//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! decomposing N pools
    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
       if (use_vertsoilc) then
          ptr2d => this%decomp_npools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%decomp_npools_vr_col(:,1,k)
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do

    if (use_vertsoilc) then
       ptr2d => this%ntrunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%ntrunc_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    if (use_nitrif_denitrif) then
       ! smin_no3_vr
       if (use_vertsoilc) then
          ptr2d => this%smin_no3_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%smin_no3_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_no3_vr'//' is required on an initialization dataset' )
       end if
    end if

    if (use_nitrif_denitrif) then
       ! smin_nh4
       if (use_vertsoilc) then
          ptr2d => this%smin_nh4_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%smin_nh4_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_nh4_vr'//' is required on an initialization dataset' )
       end if
    end if

    ! Set the integrated sminn based on sminn_vr, as is done in CNSummaryMod (this may
    ! not be the most appropriate method or place to do this)

    this%sminn_col(bounds%begc:bounds%endc) = 0._r8
    do j = 1, nlevdecomp
       do c = bounds%begc, bounds%endc
          this%sminn_col(c) = &
               this%sminn_col(c) + &
               this%sminn_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! decomp_cascade_state - the purpose of this is to check to make sure the bgc used 
    ! matches what the restart file was generated with.  
    ! add info about the SOM decomposition cascade

    if (use_century_decomp) then
       decomp_cascade_state = 1
    else
       decomp_cascade_state = 0
    end if
    ! add info about the nitrification / denitrification state
    if (use_nitrif_denitrif) then
       decomp_cascade_state = decomp_cascade_state + 10
    end if
    if (flag == 'write') itemp = decomp_cascade_state    
    call restartvar(ncid=ncid, flag=flag, varname='decomp_cascade_state', xtype=ncd_int,  &
         long_name='BGC of the model that wrote this restart file:' &
         // '  1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
         // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification', units='', &
         interpinic_flag='skip', readvar=readvar, data=itemp)
    if (flag=='read') then
       if (.not. readvar) then
          ! assume, for sake of backwards compatibility, that if decomp_cascade_state 
          ! is not in the restart file, then the current model state is the same as 
          ! the prior model state
          restart_file_decomp_cascade_state = decomp_cascade_state
          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not ' &
               // ' contain info on decomp_cascade_state used to generate the restart file.  '
          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', decomp_cascade_state
       else
          restart_file_decomp_cascade_state = itemp  
          if (decomp_cascade_state /= restart_file_decomp_cascade_state ) then
             if ( masterproc ) then
                write(iulog,*) 'CNRest: ERROR--the decomposition cascade differs between the current ' &
                     // ' model state and the model that wrote the restart file. '
                write(iulog,*) 'The model will be horribly out of equilibrium until after a lengthy spinup. '
                write(iulog,*) 'Stopping here since this is probably an error in configuring the run. '
                write(iulog,*) 'If you really wish to proceed, then override by setting '
                write(iulog,*) 'override_bgc_restart_mismatch_dump to .true. in the namelist'
                if ( .not. override_bgc_restart_mismatch_dump ) then
                   call endrun(msg= ' CNRest: Stopping. Decomposition cascade mismatch error.'//&
                        errMsg(__FILE__, __LINE__))
                endif
             endif
          endif
       end if
    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------

    ! Do nothing for write
    ! Note that the call to write spinup_state out was done in soilbiogeochem_carbonstate_inst and
    ! cannot be called again because it will try to define the variable twice
    ! when the flag below is set to define
    if (flag == 'read') then
       call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup', units='', &
            interpinic_flag='copy', readvar=readvar,  data=idata)
       if (readvar) then
          restart_file_spinup_state = idata
       else
          ! assume, for sake of backwards compatibility, that if spinup_state is not in 
          ! the restart file then current model state is the same as prior model state
          restart_file_spinup_state = spinup_state
          if ( masterproc ) then
             write(iulog,*) ' WARNING!  Restart file does not contain info ' &
                  // ' on spinup state used to generate the restart file. '
             write(iulog,*) '   Assuming the same as current setting: ', spinup_state
          end if
       end if
    end if

    ! now compare the model and restart file spinup states, and either take the 
    ! model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing pool 
    ! by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing pool 
    ! by the associated AD factor.
    ! only allow this to occur on first timestep of model run.

    if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
       if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
          if ( masterproc ) write(iulog,*) ' NitrogenStateType Restart: taking SOM pools out of AD spinup mode'
          exit_spinup = .true.
       else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
          if ( masterproc ) write(iulog,*) ' NitrogenStateType Restart: taking SOM pools into AD spinup mode'
          enter_spinup = .true.
       else
          call endrun(msg=' Error in entering/exiting spinup.  spinup_state ' &
               // ' != restart_file_spinup_state, but do not know what to do'//&
               errMsg(__FILE__, __LINE__))
       end if
       if (get_nstep() >= 2) then
          call endrun(msg=' Error in entering/exiting spinup - should occur only when nstep = 1'//&
               errMsg(__FILE__, __LINE__))
       endif
       do k = 1, ndecomp_pools
          if ( exit_spinup ) then
             m = decomp_cascade_con%spinup_factor(k)
          else if ( enter_spinup ) then
             m = 1. / decomp_cascade_con%spinup_factor(k)
          end if
          do c = bounds%begc, bounds%endc
             do j = 1, nlevdecomp
                this%decomp_npools_vr_col(c,j,k) = this%decomp_npools_vr_col(c,j,k) * m
             end do
          end do
       end do
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, num_column, filter_column, value_column )
    !
    ! !DESCRIPTION:
    ! Set nitrogen state variables
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenstate_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k      ! indices
    !------------------------------------------------------------------------

    do fi = 1,num_column
       i = filter_column(fi)

       this%sminn_col(i)       = value_column
       this%ntrunc_col(i)      = value_column
       this%cwdn_col(i)        = value_column
       if (use_nitrif_denitrif) then
          this%smin_no3_col(i) = value_column
          this%smin_nh4_col(i) = value_column
       end if
       this%totlitn_col(i)     = value_column
       this%totsomn_col(i)     = value_column
       this%totsomn_1m_col(i)  = value_column
       this%totlitn_1m_col(i)  = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%sminn_vr_col(i,j)       = value_column
          this%ntrunc_vr_col(i,j)      = value_column
          if (use_nitrif_denitrif) then
             this%smin_no3_vr_col(i,j) = value_column
             this%smin_nh4_vr_col(i,j) = value_column
          end if
       end do
    end do

    ! column and decomp_pools
    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_npools_col(i,k)    = value_column
          this%decomp_npools_1m_col(i,k) = value_column
       end do
    end do

    ! column levdecomp, and decomp_pools
    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc)
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l     ! indices
    integer  :: fc          ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

   ! vertically integrate NO3 NH4 N2O pools
   if (use_nitrif_denitrif) then
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%smin_no3_col(c) = 0._r8
         this%smin_nh4_col(c) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%smin_no3_col(c) = &
                 this%smin_no3_col(c) + &
                 this%smin_no3_vr_col(c,j) * dzsoi_decomp(j)
            
            this%smin_nh4_col(c) = &
                 this%smin_nh4_col(c) + &
                 this%smin_nh4_vr_col(c,j) * dzsoi_decomp(j)
          end do 
       end do

    end if

   ! vertically integrate each of the decomposing N pools
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%decomp_npools_col(c,l) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_npools_col(c,l) = &
                 this%decomp_npools_col(c,l) + &
                 this%decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp > 1) then

      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_npools_1m_col(c,l) = 0._r8
         end do
      end do

      ! vertically integrate each of the decomposing n pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) <= maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_npools_1m_col(c,l) = &
                       this%decomp_npools_1m_col(c,l) + &
                       this%decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) < maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_npools_1m_col(c,l) = &
                       this%decomp_npools_1m_col(c,l) + &
                       this%decomp_npools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do
      
      ! total litter nitrogen to 1 meter (TOTLITN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totlitn_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totlitn_1m_col(c) = &
                    this%totlitn_1m_col(c) + &
                    this%decomp_npools_1m_col(c,l)
            end do
         end if
      end do
      
      ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totsomn_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totsomn_1m_col(c) = this%totsomn_1m_col(c) + &
                    this%decomp_npools_1m_col(c,l)
            end do
         end if
      end do
      
   endif
   
   ! total litter nitrogen (TOTLITN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totlitn_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totlitn_col(c) = &
                 this%totlitn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do
   
   ! total soil organic matter nitrogen (TOTSOMN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totsomn_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_soil(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totsomn_col(c) = this%totsomn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do
   
   ! total cwdn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%cwdn_col(c) = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%cwdn_col(c) = this%cwdn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do

   ! total sminn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%sminn_col(c)      = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%sminn_col(c) = this%sminn_col(c) + &
              this%sminn_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! total col_ntrunc
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%ntrunc_col(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%ntrunc_col(c) = this%ntrunc_col(c) + &
              this%ntrunc_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

 end subroutine Summary

end module SoilBiogeochemNitrogenStateType
