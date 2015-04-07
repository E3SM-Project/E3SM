module SoilBiogeochemCarbonFluxType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use decompMod                          , only : bounds_type
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevgrnd, nlevdecomp
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp
  use landunit_varcon                    , only : istsoil, istcrop, istdlak 
  use ch4varcon                          , only : allowlakeprod
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use ColumnType                         , only : col                
  use LandunitType                       , only : lun                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: soilbiogeochem_carbonflux_type

     ! fire fluxes
     real(r8), pointer :: lf_conv_cflux_col                         (:)     ! (gC/m2/s) conversion C flux due to BET and BDT area decreasing (immediate loss to atm)
     real(r8), pointer :: somc_fire_col                             (:)     ! (gC/m2/s) carbon emissions due to peat burning

     ! decomposition fluxes
     real(r8), pointer :: decomp_cpools_sourcesink_col              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)  
     real(r8), pointer :: decomp_cascade_hr_vr_col                  (:,:,:) ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     real(r8), pointer :: decomp_cascade_hr_col                     (:,:)   ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
     real(r8), pointer :: decomp_cascade_ctransfer_vr_col           (:,:,:) ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_ctransfer_col              (:,:)   ! vertically-integrated (diagnostic) C transferred along decomposition cascade (gC/m2/s)
     real(r8), pointer :: decomp_k_col                              (:,:,:) ! rate constant for decomposition (1./sec)
     real(r8), pointer :: hr_vr_col                                 (:,:)   ! (gC/m3/s) total vertically-resolved het. resp. from decomposing C pools 
     real(r8), pointer :: o_scalar_col                              (:,:)   ! fraction by which decomposition is limited by anoxia
     real(r8), pointer :: w_scalar_col                              (:,:)   ! fraction by which decomposition is limited by moisture availability
     real(r8), pointer :: t_scalar_col                              (:,:)   ! fraction by which decomposition is limited by temperature
     real(r8), pointer :: som_c_leached_col                         (:)     ! (gC/m^2/s) total SOM C loss from vertical transport 
     real(r8), pointer :: decomp_cpools_leached_col                 (:,:)   ! (gC/m^2/s) C loss from vertical transport from each decomposing C pool 
     real(r8), pointer :: decomp_cpools_transport_tendency_col      (:,:,:) ! (gC/m^3/s) C tendency due to vertical transport in decomposing C pools 

     ! nitrif_denitrif
     real(r8), pointer :: phr_vr_col                                (:,:)   ! (gC/m3/s) potential hr (not N-limited) 
     real(r8), pointer :: fphr_col                                  (:,:)   ! fraction of potential heterotrophic respiration

     real(r8), pointer :: hr_col                                    (:)     ! (gC/m2/s) total heterotrophic respiration
     real(r8), pointer :: lithr_col                                 (:)     ! (gC/m2/s) litter heterotrophic respiration 
     real(r8), pointer :: somhr_col                                 (:)     ! (gC/m2/s) soil organic matter heterotrophic res   

   contains

     procedure , public  :: Init   
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory
     procedure , private :: InitCold
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary

  end type soilbiogeochem_carbonflux_type
  !------------------------------------------------------------------------

contains
   
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type)

     class(soilbiogeochem_carbonflux_type) :: this
     type(bounds_type), intent(in) :: bounds  
     character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']

     call this%InitAllocate ( bounds)
     call this%InitHistory ( bounds, carbon_type )
     call this%InitCold (bounds )

   end subroutine Init

   !------------------------------------------------------------------------
   subroutine InitAllocate(this, bounds)
     !
     ! !ARGUMENTS:
     class (soilbiogeochem_carbonflux_type) :: this 
     type(bounds_type), intent(in)    :: bounds 
     !
     ! !LOCAL VARIABLES:
     integer           :: begp,endp
     integer           :: begc,endc
     !------------------------------------------------------------------------

     begp = bounds%begp; endp = bounds%endp
     begc = bounds%begc; endc = bounds%endc

     allocate(this%t_scalar_col      (begc:endc,1:nlevdecomp_full)); this%t_scalar_col      (:,:) =spval
     allocate(this%w_scalar_col      (begc:endc,1:nlevdecomp_full)); this%w_scalar_col      (:,:) =spval
     allocate(this%o_scalar_col      (begc:endc,1:nlevdecomp_full)); this%o_scalar_col      (:,:) =spval
     allocate(this%phr_vr_col        (begc:endc,1:nlevdecomp_full)); this%phr_vr_col        (:,:) =nan 
     allocate(this%fphr_col          (begc:endc,1:nlevgrnd))       ; this%fphr_col          (:,:) =nan 
     allocate(this%som_c_leached_col (begc:endc))                  ; this%som_c_leached_col (:)   =nan
     allocate(this%somc_fire_col     (begc:endc))                  ; this%somc_fire_col     (:)   =nan
     allocate(this%lf_conv_cflux_col (begc:endc))                  ; this%lf_conv_cflux_col (:)   =nan
     allocate(this%hr_vr_col         (begc:endc,1:nlevdecomp_full)); this%hr_vr_col         (:,:) =nan

     allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
     this%decomp_cpools_sourcesink_col(:,:,:)= nan

     allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
     this%decomp_cascade_hr_vr_col(:,:,:)= spval

     allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
     this%decomp_cascade_hr_col(:,:)= nan

     allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
     this%decomp_cascade_ctransfer_vr_col(:,:,:)= nan

     allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
     this%decomp_cascade_ctransfer_col(:,:)= nan

     allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
     this%decomp_k_col(:,:,:)= spval

     allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
     this%decomp_cpools_leached_col(:,:)= nan

     allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
     this%decomp_cpools_transport_tendency_col(:,:,:)= nan

     allocate(this%hr_col                  (begc:endc)) ; this%hr_col                  (:) = nan
     allocate(this%lithr_col               (begc:endc)) ; this%lithr_col               (:) = nan
     allocate(this%somhr_col               (begc:endc)) ; this%somhr_col               (:) = nan

   end subroutine InitAllocate 

   !------------------------------------------------------------------------
   subroutine InitHistory(this, bounds, carbon_type)
     !
     ! !DESCRIPTION:
     ! add history fields for all CN variables, always set as default='inactive'
     !
     ! !USES:
     use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
     use clm_varpar , only : nlevdecomp, nlevdecomp_full
     use clm_varctl , only : hist_wrtch4diag
     use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
     !
     ! !ARGUMENTS:
     class(soilbiogeochem_carbonflux_type) :: this    
     type(bounds_type)         , intent(in) :: bounds 
     character(len=3)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
     !
     ! !LOCAL VARIABLES:
     integer           :: k,l,ii,jj,c
     character(8)      :: vr_suffix
     character(10)     :: active
     integer           :: begp,endp
     integer           :: begc,endc
     character(24)     :: fieldname
     character(100)    :: longname
     real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
     real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
     !---------------------------------------------------------------------

     begp = bounds%begp; endp = bounds%endp
     begc = bounds%begc; endc = bounds%endc

     if (nlevdecomp > 1) then
        vr_suffix = "_vr"
     else 
        vr_suffix = ""
     endif

     !-------------------------------
     ! C flux variables - native to column 
     !-------------------------------

     ! add history fields for all CLAMP CN variables

     if (carbon_type == 'c12') then

        this%hr_col(begc:endc) = spval
        call hist_addfld1d (fname='HR', units='gC/m^2/s', &
             avgflag='A', long_name='total heterotrophic respiration', &
             ptr_col=this%hr_col)

        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='LITTERC_HR', units='gC/m^2/s', &
             avgflag='A', long_name='litter C heterotrophic respiration', &
             ptr_col=this%lithr_col)
        call hist_addfld1d (fname='LITHR', units='gC/m^2/s', &
             avgflag='A', long_name='litter heterotrophic respiration', &
             ptr_col=this%lithr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='SOILC_HR', units='gC/m^2/s', &
             avgflag='A', long_name='soil C heterotrophic respiration', &
             ptr_col=this%somhr_col)
        call hist_addfld1d (fname='SOMHR', units='gC/m^2/s', &
             avgflag='A', long_name='soil organic matter heterotrophic respiration', &
             ptr_col=this%somhr_col)
        call hist_addfld1d (fname='SOILC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='soil C loss', &
             ptr_col=this%somhr_col)

        if (hist_wrtch4diag) then
           this%fphr_col(begc:endc,1:nlevgrnd) = spval
           call hist_addfld_decomp (fname='FPHR'//trim(vr_suffix), units='unitless', type2d='levdcmp', &
                avgflag='A', long_name='fraction of potential HR due to N limitation', &
                ptr_col=this%fphr_col)
        end if

        this%lf_conv_cflux_col(begc:endc) = spval
        call hist_addfld1d (fname='LF_CONV_CFLUX', units='gC/m^2/s', &
             avgflag='A', long_name='conversion carbon due to BET and BDT area decreasing', &
             ptr_col=this%lf_conv_cflux_col)   

        this%somc_fire_col(begc:endc) = spval
        call hist_addfld1d (fname='SOMC_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='C loss due to peat burning', &
             ptr_col=this%somc_fire_col)

        do k = 1, ndecomp_pools
           ! decomposition k
           data2dptr => this%decomp_k_col(:,:,k)
           fieldname = 'K_'//trim(decomp_cascade_con%decomp_pool_name_history(k))
           longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' potential loss coefficient'
           call hist_addfld_decomp (fname=fieldname, units='1/s',  type2d='levdcmp', &
                avgflag='A', long_name=longname, &
                ptr_col=data2dptr, default='inactive')
        end do

        this%decomp_cascade_hr_col(begc:endc,:)             = spval
        this%decomp_cascade_hr_vr_col(begc:endc,:,:)        = spval
        this%decomp_cascade_ctransfer_col(begc:endc,:)      = spval
        this%decomp_cascade_ctransfer_vr_col(begc:endc,:,:) = spval
        do l = 1, ndecomp_cascade_transitions

           ! output the vertically integrated fluxes only as  default
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data1dptr => this%decomp_cascade_hr_col(:,l)
              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = &
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR'
              else
                 fieldname = &
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))
              endif
              longname =  'Het. Resp. from '//&
                   trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s',  &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr)
           endif

           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data1dptr => this%decomp_cascade_ctransfer_col(:,l)
              fieldname = &
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'C'
              longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                   ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr)
           endif

           ! output the vertically resolved fluxes 
           if ( nlevdecomp_full > 1 ) then  
              !-- HR fluxes (none from CWD)
              if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
                 data2dptr => this%decomp_cascade_hr_vr_col(:,:,l)
                 ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
                 ii = 0
                 do jj = 1, ndecomp_cascade_transitions
                    if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
                 end do
                 if ( ii == 1 ) then
                    fieldname = &
                         trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                         //'_HR'//trim(vr_suffix)
                 else
                    fieldname = &
                         trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                         trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))&
                         //trim(vr_suffix)
                 endif
                 longname =  'Het. Resp. from '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
                 call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              endif

              !-- transfer fluxes (none from terminal pool, if present)
              if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                 data2dptr => this%decomp_cascade_ctransfer_vr_col(:,:,l)
                 fieldname = &
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                      //'C'//trim(vr_suffix)
                 longname =  'decomp. of '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                      ' C to '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
                 call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              endif
           end if

        end do

        this%t_scalar_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='T_SCALAR', units='unitless',  type2d='levdcmp', &
             avgflag='A', long_name='temperature inhibition of decomposition', &
             ptr_col=this%t_scalar_col)

        this%w_scalar_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='W_SCALAR', units='unitless',  type2d='levdcmp', &
             avgflag='A', long_name='Moisture (dryness) inhibition of decomposition', &
             ptr_col=this%w_scalar_col)

        this%o_scalar_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='O_SCALAR', units='unitless', type2d='levdcmp', &
             avgflag='A', long_name='fraction by which decomposition is reduced due to anoxia', &
             ptr_col=this%o_scalar_col)

        this%som_c_leached_col(begc:endc) = spval
        call hist_addfld1d (fname='SOM_C_LEACHED', units='gC/m^2/s', &
             avgflag='A', long_name='total flux of C from SOM pools due to leaching', &
             ptr_col=this%som_c_leached_col)!, default='inactive')

        this%decomp_cpools_leached_col(begc:endc,:) = spval
        this%decomp_cpools_transport_tendency_col(begc:endc,:,:) = spval
        do k = 1, ndecomp_pools
           if ( .not. decomp_cascade_con%is_cwd(k) ) then
              data1dptr => this%decomp_cpools_leached_col(:,k)
              fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_LEACHING'
              longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C leaching loss'
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr)!, default='inactive')

              data2dptr => this%decomp_cpools_transport_tendency_col(:,:,k)
              fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TNDNCY_VERT_TRANSPORT'
              longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C tendency due to vertical transport'
              call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
        end do

        if ( nlevdecomp_full > 1 ) then
           this%hr_vr_col(begc:endc,:) = spval
           call hist_addfld2d (fname='HR_vr', units='gC/m^3/s', type2d='levdcmp', &
                avgflag='A', long_name='total vertically resolved heterotrophic respiration', &
                ptr_col=this%hr_vr_col)
        endif

     end if

     !-------------------------------
     ! C13 flux variables - native to column 
     !-------------------------------

     if ( carbon_type == 'c13' ) then

        this%hr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_HR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total heterotrophic respiration', &
             ptr_col=this%hr_col)

        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_LITHR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C litterfall to litter 3 C', &
             ptr_col=this%lithr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_SOMHR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 soil organic matter heterotrophic respiration', &
             ptr_col=this%somhr_col)


        this%decomp_cascade_hr_col(begc:endc,:)             = spval
        this%decomp_cascade_hr_vr_col(begc:endc,:,:)        = spval
        this%decomp_cascade_ctransfer_col(begc:endc,:)      = spval
        this%decomp_cascade_ctransfer_vr_col(begc:endc,:,:) = spval
        do l = 1, ndecomp_cascade_transitions
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data2dptr => this%decomp_cascade_hr_vr_col(:,:,l)
              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR'//trim(vr_suffix)
              else
                 fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))//&
                      trim(vr_suffix)
              endif
              longname =  'C13 Het. Resp. from '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data2dptr => this%decomp_cascade_ctransfer_vr_col(:,:,l)
              fieldname = 'C13_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                   //'C_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                   //'C'//trim(vr_suffix)
              longname =  'C13 decomp. of '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))&
                   //' C to '//&
                   trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
              call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
        end do

     end if

     !-------------------------------
     ! C14 flux variables - native to column 
     !-------------------------------

     if (carbon_type == 'c14') then

        this%hr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_HR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total heterotrophic respiration', &
             ptr_col=this%hr_col)

        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_LITHR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C litterfall to litter 3 C', &
             ptr_col=this%lithr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_SOMHR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 soil organic matter heterotrophic respiration', &
             ptr_col=this%somhr_col)

        this%decomp_cascade_hr_col(begc:endc,:)             = spval
        this%decomp_cascade_hr_vr_col(begc:endc,:,:)        = spval
        this%decomp_cascade_ctransfer_col(begc:endc,:)      = spval
        this%decomp_cascade_ctransfer_vr_col(begc:endc,:,:) = spval

        do l = 1, ndecomp_cascade_transitions
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data2dptr => this%decomp_cascade_hr_vr_col(:,:,l)

              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR'//trim(vr_suffix)
              else
                 fieldname = 'C14_'//&
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))&
                      //trim(vr_suffix)
              endif
              longname =  'C14 Het. Resp. from '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif

           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data2dptr => this%decomp_cascade_ctransfer_vr_col(:,:,l)

              fieldname = 'C14_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                   //'C_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                   //'C'//trim(vr_suffix)
              longname =  'C14 decomp. of '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                   ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
              call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
        end do

     end if

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       this%fphr_col(c,nlevdecomp+1:nlevgrnd) = 0._r8 !used to be in ch4Mod
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%fphr_col(c,nlevdecomp+1:nlevgrnd) = 0._r8 
       else if (lun%itype(l) == istdlak .and. allowlakeprod) then
          this%fphr_col(c,:) = spval
       else  ! Inactive CH4 columns
          this%fphr_col(c,:) = spval
       end if

       ! also initialize dynamic landcover fluxes so that they have
       ! real values on first timestep, prior to calling pftdyn_cnbal
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%lf_conv_cflux_col(c) = 0._r8
       end if
    end do

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: c,l
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! initialize fields for special filters

    call this%SetValues (num_column=num_special_col, filter_column=special_col, &
         value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonflux_type) :: this
    type(bounds_type) , intent(in)        :: bounds  
    type(file_desc_t) , intent(inout)     :: ncid   ! netcdf id
    character(len=*)  , intent(in)        :: flag   !'read', 'write', 'define'
    !-----------------------------------------------------------------------

    ! Nothing for now

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon fluxes
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_carbonflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k,l    ! indices
    !------------------------------------------------------------------------

    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_hr_col(i,l)             = value_column
             this%decomp_cascade_hr_vr_col(i,j,l)        = value_column
             this%decomp_cascade_ctransfer_col(i,l)      = value_column
             this%decomp_cascade_ctransfer_vr_col(i,j,l) = value_column
             this%decomp_k_col(i,j,l)                    = value_column
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cpools_leached_col(i,k) = value_column
       end do
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_transport_tendency_col(i,j,k) = value_column
             this%decomp_cpools_sourcesink_col(i,j,k)         = value_column  
          end do
       end do
    end do

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%hr_vr_col(i,j) = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)
       this%hr_col(i)            = value_column
       this%somc_fire_col(i)     = value_column  
       this%som_c_leached_col(i) = value_column
       this%somhr_col(i)         = value_column
       this%lithr_col(i)         = value_column
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonflux_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j          ! indices
    !-----------------------------------------------------------------------

    ! set column-level conversion and product pool fluxes
    ! to 0 at the beginning of every timestep

    do c = bounds%begc,bounds%endc
       this%lf_conv_cflux_col(c) = 0._r8
    end do


  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonflux_type)           :: this
    type(bounds_type)               , intent(in)    :: bounds          
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l
    integer  :: fc       
    !-----------------------------------------------------------------------

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%som_c_leached_col(c) = 0._r8
    end do

    ! vertically integrate HR and decomposition cascade fluxes
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cascade_hr_col(c,k) = &
                  this%decomp_cascade_hr_col(c,k) + &
                  this%decomp_cascade_hr_vr_col(c,j,k) * dzsoi_decomp(j) 

             this%decomp_cascade_ctransfer_col(c,k) = &
                  this%decomp_cascade_ctransfer_col(c,k) + &
                  this%decomp_cascade_ctransfer_vr_col(c,j,k) * dzsoi_decomp(j) 
          end do
       end do
    end do

    ! total heterotrophic respiration, vertically resolved (HR)
    do j = 1,nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%hr_vr_col(c,j) = 0._r8
       end do
    end do
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%hr_vr_col(c,j) = &
                  this%hr_vr_col(c,j) + &
                  this%decomp_cascade_hr_vr_col(c,j,k)
          end do
       end do
    end do

    ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_cpools_leached_col(c,l) = 0._r8
       end do
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_leached_col(c,l) = this%decomp_cpools_leached_col(c,l) + &
                  this%decomp_cpools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%som_c_leached_col(c) = this%som_c_leached_col(c) + this%decomp_cpools_leached_col(c,l)
       end do
    end do

    ! soil organic matter heterotrophic respiration 
    associate(is_soil => decomp_cascade_con%is_soil) ! TRUE => pool is a soil pool  
      do k = 1, ndecomp_cascade_transitions
         if ( is_soil(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%somhr_col(c) = this%somhr_col(c) + this%decomp_cascade_hr_col(c,k)
            end do
         end if
      end do
    end associate

    ! litter heterotrophic respiration (LITHR)
    associate(is_litter => decomp_cascade_con%is_litter) ! TRUE => pool is a litter pool
      do k = 1, ndecomp_cascade_transitions
         if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%lithr_col(c) = this%lithr_col(c) + this%decomp_cascade_hr_col(c,k)
            end do
         end if
      end do
    end associate

    ! total heterotrophic respiration (HR)
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%hr_col(c) = &
            this%lithr_col(c) + &
            this%somhr_col(c)
    end do

  end subroutine Summary

end module SoilBiogeochemCarbonFluxType


