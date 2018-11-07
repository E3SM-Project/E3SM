module CNBalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon mass balance checking.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog, use_nitrif_denitrif, use_fates
  use clm_time_manager    , only : get_step_size,get_nstep
  use clm_varpar          , only : crop_prog
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use ColumnType          , only : col_pp                
  use GridcellType        , only : grc_pp
  use clm_varpar          , only : nlevdecomp
  use clm_varcon          , only : dzsoi_decomp
  use clm_varctl          , only : nu_com
  use clm_varctl          , only : ECA_Pconst_RGspin

  use spmdMod             , only : iam
  use CNDecompCascadeConType , only : decomp_cascade_con
  use clm_varpar          , only: ndecomp_cascade_transitions
  use subgridAveMod       , only : p2c
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode, pf_hmode
  ! forest fertilization experiment
  use clm_time_manager    , only : get_curr_date
  use CNStateType         , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use clm_varctl          , only : forest_fert_exp
  use VegetationType      , only : veg_pp
  use pftvarcon           , only: noveg
  use clm_varctl          , only : NFIX_PTASE_plant

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginCBalance
  public :: BeginNBalance
  public :: BeginPBalance
  public :: CBalanceCheck
  public :: NBalanceCheck
  public :: PBalanceCheck
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BeginCBalance(bounds, num_soilc, filter_soilc, &
       carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning carbon balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                        & 
         totcolc   =>  carbonstate_vars%totcolc_col , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         col_begcb =>  carbonstate_vars%begcb_col     & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
         )

      ! calculate beginning column-level carbon balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begcb(c) = totcolc(c)
!         if(c==5657 .and. .false.)then
!           print*,'totbgc blgc',carbonstate_vars%totabgc_col(c),carbonstate_vars%totblgc_col(c)
         carbonstate_vars%begabgc_col(c) = carbonstate_vars%totabgc_col(c)
         carbonstate_vars%begblgc_col(c) = carbonstate_vars%totblgc_col(c)
!         endif
!         col_pp%debug_flag(c)=(get_nstep()>=412) .and. (c==20024)
!         if(c==20024)print*,'flag',col_pp%debug_flag(c),get_nstep()
      end do
!      if(iam==6)col_pp%debug_flag(17236)=.true.
    end associate

  end subroutine BeginCBalance
 
  !-----------------------------------------------------------------------
  subroutine BeginNBalance(bounds, num_soilc, filter_soilc, &
       nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning nitrogen balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds          
    integer                  , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                         & 
         totcoln   => nitrogenstate_vars%totcoln_col , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         col_begnb => nitrogenstate_vars%begnb_col     & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         )

      ! calculate beginning column-level nitrogen balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begnb(c) = totcoln(c)
         nitrogenstate_vars%begabgn_col(c)=nitrogenstate_vars%totabgn_col(c)
         nitrogenstate_vars%begblgn_col(c)=nitrogenstate_vars%totblgn_col(c)
!         col_pp%debug_flag(c) = (abs(grc_pp%latdeg(col_pp%gridcell(c))-82.4210526315788)+abs(grc_pp%londeg(col_pp%gridcell(c))-317.50) < 1.e-5 &
!            .and. get_nstep()==2)
      end do

    end associate

  end subroutine BeginNBalance

  !-----------------------------------------------------------------------
  subroutine BeginPBalance(bounds, num_soilc, filter_soilc, &
       phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning phosphorus balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds          
    integer                  , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                           &
         totcolp   => phosphorusstate_vars%totcolp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         !X.YANG - checking P balance problem, starting from VEGP 
         totpftp   => phosphorusstate_vars%totpftp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp   => phosphorusstate_vars%totsomp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp   => phosphorusstate_vars%cwdp_col       , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp   => phosphorusstate_vars%totlitp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         sminp   => phosphorusstate_vars%sminp_col     , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
 
         col_begpb => phosphorusstate_vars%begpb_col     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         )

      ! calculate beginning column-level phosphorus balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begpb(c) = totcolp(c)
         phosphorusstate_vars%begabgp_col(c) = phosphorusstate_vars%totabgp_col(c)
         phosphorusstate_vars%begblgp_col(c) = phosphorusstate_vars%totblgp_col(c)
!         col_begpb(c) = totpftp(c)
!         col_begpb(c) = sminp(c)
!         col_begpb(c) = totsomp(c)+totlitp(c)+cwdp(c)
!        if(c==1)then
!          print*,'============================='
!          print*,'nstep,tbegp=',c,get_nstep()
!          print*, 'occl p=', phosphorusstate_vars%occlp_col(c)
!          print*,'begp   =', col_begpb(c)
!          print*,'begabgp=', phosphorusstate_vars%totabgp_col(c),phosphorusstate_vars%occlp_col(c)
!          print*,'begblgp=', phosphorusstate_vars%totblgp_col(c)+phosphorusstate_vars%occlp_col(c)
!        endif
      end do

    end associate

  end subroutine BeginPBalance

  !-----------------------------------------------------------------------
  subroutine CBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       carbonstate_vars, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform carbon mass conservation check for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,err_index    ! indices
    integer  :: fc             ! lake filter indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: col_cinputs
    real(r8) :: col_coutputs
    real(r8) :: err_blg,err_abg
    !-----------------------------------------------------------------------

    associate(                                                                   & 
         totcolc                 =>    carbonstate_vars%totcolc_col            , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         
         gpp                     =>    carbonflux_vars%gpp_col                 , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production      
         er                      =>    carbonflux_vars%er_col                  , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs          =>    carbonflux_vars%fire_closs_col          , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_hrv_xsmrpool_to_atm =>    carbonflux_vars%hrv_xsmrpool_to_atm_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality 
         dwt_closs               =>    carbonflux_vars%dwt_closs_col           , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total carbon loss from product pools and conversion
         product_closs           =>    carbonflux_vars%product_closs_col       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total wood product carbon loss
         som_c_leached           =>    carbonflux_vars%som_c_leached_col       , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss from vertical transport 
         col_decompc_delta       =>    carbonflux_vars%externalc_to_decomp_delta_col, & ! Input: [real(r8) (:) ] (gC/m2/s) summarized net change of whole column C i/o to decomposing pool bwtn time-step
         
         col_begcb               =>    carbonstate_vars%begcb_col              , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         col_endcb               =>    carbonstate_vars%endcb_col              , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2) 
         col_errcb               =>    carbonstate_vars%errcb_col              , & ! Output: [real(r8) (:) ]  carbon balance error for the timestep (gC/m**2)
         som_c_runoff_col        =>    carbonflux_vars%som_c_runoff_col          &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! calculate the total column-level carbon storage, for mass conservation check
         col_endcb(c) = totcolc(c)

         ! calculate total column-level inputs
         col_cinputs = gpp(c)

         ! calculate total column-level outputs
         ! er = ar + hr, col_fire_closs includes pft-level fire losses
         col_coutputs = er(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c) + col_hrv_xsmrpool_to_atm(c)

         ! subtract leaching flux
         col_coutputs = col_coutputs - som_c_leached(c) + som_c_runoff_col(c)

         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs - col_coutputs)*dt - (col_endcb(c) - col_begcb(c))

         ! adjusting the time-lag of org. C increments to decomposing pools when coupled with PFLOTRAN bgc
         ! (because PF bgc uses the extern C as sink (in, + ) at previous time-step,
         ! but note that it includes possible negative adding)
         if (use_pflotran .and. pf_cmode) then
            col_errcb(c) = col_errcb(c) - col_decompc_delta(c)*dt
            ! here is '-' adjustment. It says that the adding to PF decomp c pools was less.
         end if

         ! check for significant errors
         if (abs(col_errcb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
            exit
         end if
      end do ! end of columns loop
      
      ! Consider adapting this check to be fates compliant (rgk 04-2017)
      if (.not. use_fates) then
         if (err_found) then
            c = err_index
            write(iulog,*)'column cbalance error = ', col_errcb(c), c, get_nstep()
            write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
            write(iulog,*)'input                 = ',col_cinputs*dt
            write(iulog,*)'output                = ',col_coutputs*dt
            write(iulog,*)'er                    = ',er(c)*dt
            write(iulog,*)'hr                    = ',carbonflux_vars%hr_col(c)*dt
            write(iulog,*)'fire                  = ',col_fire_closs(c)*dt
            write(iulog,*)'dwt                   = ',dwt_closs(c)*dt
            write(iulog,*)'product               = ',product_closs(c)*dt
            write(iulog,*)'hrv                   = ',col_hrv_xsmrpool_to_atm(c)*dt
            write(iulog,*)'leach                 = ',som_c_leached(c)*dt
            write(iulog,*)'runoff som c          = ',som_c_runoff_col(c)*dt
            write(iulog,*)'begcb                 = ',col_begcb(c)
            write(iulog,*)'endcb                 = ',col_endcb(c)
            write(iulog,*)'begabg                = ',carbonstate_vars%begabgc_col(c)
            write(iulog,*)'begblg                = ',carbonstate_vars%begblgc_col(c)
            write(iulog,*)'delta store           = ',col_endcb(c)-col_begcb(c)
            write(iulog,*)'abgc                  = ',carbonstate_vars%totabgc_col(c)
            write(iulog,*)'blgc                  = ',carbonstate_vars%totblgc_col(c)
            write(iulog,*)'abg_fire_loss         = ',carbonflux_vars%fire_closs_p2c_col(c)*dt
            write(iulog,*)'blg_fire_loss         = ',carbonflux_vars%fire_decomp_closs_col(c)*dt
            write(iulog,*)'ar                    = ',carbonflux_vars%ar_col(c)*dt
            write(iulog,*)'plt2soi               = ',carbonflux_vars%cflx_plant_to_soilbgc_col(c)*dt
            if (use_pflotran .and. pf_cmode) then
               write(iulog,*)'pf_delta_decompc      = ',col_decompc_delta(c)*dt
            end if
            err_blg=carbonstate_vars%begblgc_col(c)+carbonflux_vars%cflx_plant_to_soilbgc_col(c)*dt-&
               carbonflux_vars%fire_decomp_closs_col(c)*dt- &
               carbonflux_vars%hr_col(c)*dt-som_c_leached(c)*dt-som_c_runoff_col(c)*dt-carbonstate_vars%totblgc_col(c)
            write(iulog,*)'delta_blg             = ',carbonstate_vars%begblgc_col(c)-carbonstate_vars%totblgc_col(c)
            write(iulog,*)'errblg                = ',err_blg
            err_abg=carbonstate_vars%begabgc_col(c)-carbonstate_vars%totabgc_col(c) &
               +gpp(c)*dt-carbonflux_vars%ar_col(c)*dt &
               -carbonflux_vars%cflx_plant_to_soilbgc_col(c)*dt-dwt_closs(c)*dt-product_closs(c)*dt-col_hrv_xsmrpool_to_atm(c)*dt &
               -carbonflux_vars%fire_closs_p2c_col(c)*dt
            write(iulog,*)'delta_abg             = ',carbonstate_vars%begabgc_col(c)-carbonstate_vars%totabgc_col(c)
            write(iulog,*)'errabg                = ',err_abg
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
      end if !use_fates

    end associate

  end subroutine CBalanceCheck

  !-----------------------------------------------------------------------
  subroutine NBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform nitrogen mass conservation check
    ! for column and pft
    !
    use tracer_varcon,  only : is_active_betr_bgc
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds          
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j,p  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)

    integer:: kyr                     ! current year 
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31) 
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day) 
    real(r8):: err_blg, err_abg
    !-----------------------------------------------------------------------

    associate(                                                                 &
         totcoln             =>    nitrogenstate_vars%totcoln_col            , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         ndep_to_sminn       =>    nitrogenflux_vars%ndep_to_sminn_col       , & ! Input:  [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         nfix_to_sminn       =>    nitrogenflux_vars%nfix_to_sminn_col       , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         nfix_to_ecosysn     =>    nitrogenflux_vars%nfix_to_ecosysn_col       , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         fert_to_sminn       =>    nitrogenflux_vars%fert_to_sminn_col       , & ! Input:  [real(r8) (:)]                                          
         soyfixn_to_sminn    =>    nitrogenflux_vars%soyfixn_to_sminn_col    , & ! Input:  [real(r8) (:)]                                          
         supplement_to_sminn =>    nitrogenflux_vars%supplement_to_sminn_col , & ! Input:  [real(r8) (:)]  supplemental N supply (gN/m2/s)         
         denit               =>    nitrogenflux_vars%denit_col               , & ! Input:  [real(r8) (:)]  total rate of denitrification (gN/m2/s) 
         sminn_leached       =>    nitrogenflux_vars%sminn_leached_col       , & ! Input:  [real(r8) (:)]  soil mineral N pool loss to leaching (gN/m2/s)
         smin_no3_leached    =>    nitrogenflux_vars%smin_no3_leached_col    , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to leaching (gN/m2/s)
         smin_no3_runoff     =>    nitrogenflux_vars%smin_no3_runoff_col     , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to runoff (gN/m2/s)
         smin_nh4_leached    =>    nitrogenflux_vars%smin_nh4_leached_col    , & ! Input:  [real(r8) (:)]  soil mineral nh4 pool loss to leaching (gN/m2/s)
         smin_nh4_runoff     =>    nitrogenflux_vars%smin_nh4_runoff_col     , & ! Input:  [real(r8) (:)]  soil mineral nh4 pool loss to runoff (gN/m2/s)
         f_n2o_nit           =>    nitrogenflux_vars%f_n2o_nit_col           , & ! Input:  [real(r8) (:)]  flux of N2o from nitrification [gN/m^2/s]
         col_fire_nloss      =>    nitrogenflux_vars%fire_nloss_col          , & ! Input:  [real(r8) (:)]  total column-level fire N loss (gN/m2/s)
         dwt_nloss           =>    nitrogenflux_vars%dwt_nloss_col           , & ! Input:  [real(r8) (:)]  (gN/m2/s) total nitrogen loss from product pools and conversion
         product_nloss       =>    nitrogenflux_vars%product_nloss_col       , & ! Input:  [real(r8) (:)]  (gN/m2/s) total wood product nitrogen loss
         som_n_leached       =>    nitrogenflux_vars%som_n_leached_col       , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport
         supplement_to_plantn=> nitrogenflux_vars%supplement_to_plantn       , &
         ! pflotran:
         col_decompn_delta   =>    nitrogenflux_vars%externaln_to_decomp_delta_col  , & ! Input: [real(r8) (:) ] (gN/m2/s) summarized net change of whole column N i/o to decomposing pool bwtn time-step

         col_ninputs         =>    nitrogenflux_vars%ninputs_col             , & ! Output: [real(r8) (:)]  column-level N inputs (gN/m2/s)
         col_noutputs        =>    nitrogenflux_vars%noutputs_col            , & ! Output: [real(r8) (:)]  column-level N outputs (gN/m2/s)
         col_begnb           =>    nitrogenstate_vars%begnb_col              , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         col_endnb           =>    nitrogenstate_vars%endnb_col              , & ! Output: [real(r8) (:)]  nitrogen mass, end of time step (gN/m**2)
         col_errnb           =>    nitrogenstate_vars%errnb_col              , & ! Output: [real(r8) (:)]  nitrogen balance error for the timestep (gN/m**2)
         som_n_runoff_col    =>    nitrogenflux_vars%som_n_runoff_col        , &
         nh3_soi_flx_col     =>    nitrogenflux_vars%nh3_soi_flx_col           &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      call get_curr_date(kyr, kmo, kda, mcsec)

      err_found = .false.
      ! column loop
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level nitrogen storage, for mass conservation check
         col_endnb(c) = totcoln(c)

         ! calculate total column-level inputs
         if (NFIX_PTASE_plant) then
            col_ninputs(c) = ndep_to_sminn(c) + nfix_to_ecosysn(c) + supplement_to_sminn(c)
         else
            col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)
         end if

         if (crop_prog) col_ninputs(c) = col_ninputs(c) + &
              fert_to_sminn(c) + soyfixn_to_sminn(c)

         do p = col_pp%pfti(c), col_pp%pftf(c)
            if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
                col_ninputs(c) = col_ninputs(c) + supplement_to_plantn(p) * veg_pp%wtcol(p)
            end if
         end do

         ! forest fertilization
         if (forest_fert_exp) then
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 1 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               col_ninputs(c) = col_ninputs(c) + fert_dose(c,kmo)/dt
             end if
         end if

         ! calculate total column-level outputs
         col_noutputs(c) = denit(c) + col_fire_nloss(c) + dwt_nloss(c) + product_nloss(c)

         if(is_active_betr_bgc)then
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c) &
               + smin_nh4_leached(c) + smin_nh4_runoff(c)+ nh3_soi_flx_col(c)
         else
           if (.not. use_nitrif_denitrif) then
            col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
           else
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            if(use_pflotran .and. pf_cmode) then
               ! inclusion of aq. NH4 transport by PFLOTRAN-bgc
               col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
            else

               col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)

            endif

           end if
         endif

         col_noutputs(c) = col_noutputs(c) - som_n_leached(c) + som_n_runoff_col(c) 

         ! calculate the total column-level nitrogen balance error for this time step
         col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
              (col_endnb(c) - col_begnb(c))

         ! adjusting the time-lag of org. N increments to decomposing pools when coupled with PFLOTRAN bgc
         ! (because PF bgc uses the extern N sink (in, +) at previous time-step,
         ! but note that it includes possible negative adding)
         if (use_pflotran .and. pf_cmode) then
            col_errnb(c) = col_errnb(c) - col_decompn_delta(c)*dt
            ! here is '-' adjustment. It says that the adding to PF decomp n pools was less.
         end if

!         if(c==38749 .and. get_nstep()==2)then
!            err_found = .true.
!            err_index = c
!            exit
!         endif
         if (abs(col_errnb(c)) > 1e-8_r8)then ! .or. nitrogenstate_vars%sminn_col(c)>1000._r8) then
            err_found = .true.
            err_index = c
            exit
         end if


      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column nbalance error = ', col_errnb(c), c, get_nstep()
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
         write(iulog,*)'begnb                 = ',col_begnb(c)
         write(iulog,*)'endnb                 = ',col_endnb(c)
         write(iulog,*)'begabgn               = ',nitrogenstate_vars%begabgn_col(c)
         write(iulog,*)'begblgn               = ',nitrogenstate_vars%begblgn_col(c)
         write(iulog,*)'delta_store           = ',col_endnb(c)-col_begnb(c)
         write(iulog,*)'net_flux              = ',(col_ninputs(c)-col_noutputs(c))*dt
         write(iulog,*)'input_mass            = ',col_ninputs(c)*dt
         write(iulog,*)'ndep                  = ',ndep_to_sminn(c)*dt
         write(iulog,*)'nfix_to_ecosys        = ',nfix_to_ecosysn(c)*dt
         write(iulog,*)'nfix_to_soil          = ',nfix_to_sminn(c)*dt
         write(iulog,*)'nsup                  = ',supplement_to_sminn(c)*dt
         if(crop_prog) then
            write(iulog,*)'fertm                 = ',fert_to_sminn(c)*dt
            write(iulog,*)'soyfx                 = ',soyfixn_to_sminn(c)*dt
         endif
         write(iulog,*)'fire                  = ',col_fire_nloss(c)*dt
         write(iulog,*)'abg fire              = ',nitrogenflux_vars%fire_nloss_p2c_col(c)*dt
         write(iulog,*)'dwt                   = ',dwt_nloss(c)*dt
         write(iulog,*)'prod                  = ',product_nloss(c)*dt
         write(iulog,*)'output_mass           = ',col_noutputs(c)*dt
         write(iulog,*)'denit                 = ',denit(c)*dt
         write(iulog,*)'n2onit                = ',f_n2o_nit(c)*dt
         write(iulog,*)'no3_leach             = ', smin_no3_leached(c)*dt
         write(iulog,*)'no3_runof             = ', smin_no3_runoff(c)*dt
         write(iulog,*)'nh4_leach             = ',smin_nh4_leached(c)*dt
         write(iulog,*)'nh4_runof             = ',smin_nh4_runoff(c)*dt
         write(iulog,*)'som_n_leach           = ', -som_n_leached(c)*dt
         write(iulog,*)'som_n_runoff          = ',som_n_runoff_col(c)*dt
         write(iulog,*)'nh3_soile             = ',nh3_soi_flx_col(c)*dt
         write(iulog,*)'decompfireloss        = ',nitrogenflux_vars%fire_decomp_nloss_col(c)*dt
         write(iulog,*)'plt2soil              = ',nitrogenflux_vars%nflx_plant_to_soilbgc_col(c)*dt
         write(iulog,*)'soil2plt              = ',nitrogenflux_vars%sminn_to_plant_col(c)*dt
         write(iulog,*)'nsup_to_surf          = ',nitrogenflux_vars%supplement_to_sminn_surf_col(c)*dt
         write(iulog,*)'abgn                  = ',nitrogenstate_vars%totabgn_col(c)
         write(iulog,*)'blgn                  = ',nitrogenstate_vars%totblgn_col(c)
         if (use_pflotran .and. pf_cmode) then
            write(iulog,*)'pf_delta_decompn      = ',col_decompn_delta(c)*dt
         end if
         err_blg=nitrogenstate_vars%begblgn_col(c)-nitrogenstate_vars%totblgn_col(c) &
            + ndep_to_sminn(c)*dt+nfix_to_sminn(c)*dt+supplement_to_sminn(c)*dt       &
            + fert_to_sminn(c)*dt+soyfixn_to_sminn(c)*dt &
            + nitrogenflux_vars%nflx_plant_to_soilbgc_col(c)*dt &
            - denit(c)*dt-f_n2o_nit(c)*dt-smin_no3_leached(c)*dt &
            - smin_no3_runoff(c)*dt+som_n_leached(c)*dt-som_n_runoff_col(c)*dt &
            - nitrogenflux_vars%fire_decomp_nloss_col(c)*dt &
            - smin_nh4_runoff(c)*dt-smin_nh4_leached(c)*dt &
            - nitrogenflux_vars%sminn_to_plant_col(c)*dt-nh3_soi_flx_col(c)*dt
         write(iulog,*)'errblg =',err_blg
         err_abg = nitrogenstate_vars%begabgn_col(c) - nitrogenstate_vars%totabgn_col(c)  &
            + (nfix_to_ecosysn(c)-nfix_to_sminn(c))*dt &
            + nitrogenflux_vars%sminn_to_plant_col(c)*dt  &
            + nitrogenflux_vars%supplement_to_sminn_surf_col(c)*dt &
            - nitrogenflux_vars%nflx_plant_to_soilbgc_col(c)*dt &
            - dwt_nloss(c)*dt &
            - product_nloss(c)*dt &
            - nitrogenflux_vars%fire_nloss_p2c_col(c)*dt
         write(iulog,*)'errabg =',err_abg


         call endrun(msg=errMsg(__FILE__, __LINE__))

      end if
    end associate

  end subroutine NBalanceCheck


  !-----------------------------------------------------------------------
  subroutine PBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform phosphorus mass conservation check
    ! for column and pft
    !
    use CNAllocationMod         , only: suplphos
    use tracer_varcon,  only : is_active_betr_bgc
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds          
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j,k,p  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)

    real(r8) :: leafp_to_litter_col(bounds%begc:bounds%endc) 
    real(r8) :: frootp_to_litter_col(bounds%begc:bounds%endc) 
    real(r8):: flux_mineralization_col(bounds%begc:bounds%endc)   !  local temperary variable

    integer:: kyr                     ! current year 
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31) 
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day) 
    real(r8):: err_blg, err_abg
    !-----------------------------------------------------------------------

    associate(                                                                   &
         totcolp             =>    phosphorusstate_vars%totcolp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         supplement_to_sminp =>    phosphorusflux_vars%supplement_to_sminp_col , & ! Input:  [real(r8) (:)]  supplemental P supply (gP/m2/s)         
         sminp_leached       =>    phosphorusflux_vars%sminp_leached_col       , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         col_fire_ploss      =>    phosphorusflux_vars%fire_ploss_col          , & ! Input:  [real(r8) (:)]  total column-level fire P loss (gP/m2/s)
         dwt_ploss           =>    phosphorusflux_vars%dwt_ploss_col           , & ! Input:  [real(r8) (:)]  (gP/m2/s) total phosphorus loss from product pools and conversion
         product_ploss       =>    phosphorusflux_vars%product_ploss_col       , & ! Input:  [real(r8) (:)]  (gP/m2/s) total wood product phosphorus loss
         primp_to_labilep    =>    phosphorusflux_vars%primp_to_labilep_col    , &
         secondp_to_occlp    =>    phosphorusflux_vars%secondp_to_occlp_col    , &
         fert_p_to_sminp     =>    phosphorusflux_vars%fert_p_to_sminp_col     , &
         supplement_to_plantp=>    phosphorusflux_vars%supplement_to_plantp    , &
 
         col_pinputs         =>    phosphorusflux_vars%pinputs_col             , & ! Output: [real(r8) (:)]  column-level P inputs (gP/m2/s)
         col_poutputs        =>    phosphorusflux_vars%poutputs_col            , & ! Output: [real(r8) (:)]  column-level P outputs (gP/m2/s)
         col_begpb           =>    phosphorusstate_vars%begpb_col              , & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         col_endpb           =>    phosphorusstate_vars%endpb_col              , & ! Output: [real(r8) (:)]  phosphorus mass, end of time step (gP/m**2)
         col_errpb           =>    phosphorusstate_vars%errpb_col              , & ! Output: [real(r8) (:)]  phosphorus balance error for the timestep (gP/m**2)

         !X.YANG testing P Balance, from VEGP
         totpftp             =>    phosphorusstate_vars%totpftp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp             =>    phosphorusstate_vars%totsomp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp                =>    phosphorusstate_vars%cwdp_col               , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp             =>    phosphorusstate_vars%totlitp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         sminp                 =>  phosphorusstate_vars%sminp_col              , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         leafp_to_litter       =>  phosphorusflux_vars%leafp_to_litter_patch   , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         frootp_to_litter      =>  phosphorusflux_vars%frootp_to_litter_patch  , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         sminp_to_plant        =>  phosphorusflux_vars%sminp_to_plant_col      , &
         cascade_receiver_pool =>  decomp_cascade_con%cascade_receiver_pool    , &
         pf                    =>  phosphorusflux_vars                         , &
         ps                    =>  phosphorusstate_vars                        , &
         som_p_leached         =>  phosphorusflux_vars%som_p_leached_col       , &
         supplement_to_sminp_surf_col => phosphorusflux_vars%supplement_to_sminp_surf_col, &
         sminp_runoff_col      => phosphorusflux_vars%sminp_runoff_col   , &
         som_p_runoff_col      => phosphorusflux_vars%som_p_runoff_col     &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      call get_curr_date(kyr, kmo, kda, mcsec)

      err_found = .false.

      call p2c(bounds,num_soilc,filter_soilc, &
           leafp_to_litter(bounds%begp:bounds%endp), &
           leafp_to_litter_col(bounds%begc:bounds%endc))
      call p2c(bounds,num_soilc,filter_soilc, &
           frootp_to_litter(bounds%begp:bounds%endp), &
           frootp_to_litter_col(bounds%begc:bounds%endc))

      !! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes
      !! - X.YANG
        ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         flux_mineralization_col(c) = 0._r8
      enddo

      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  flux_mineralization_col(c) = flux_mineralization_col(c) - &
                                               pf%decomp_cascade_sminp_flux_col(c,k)
               end do
         else
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                               pf%decomp_cascade_sminp_flux_col(c,k)

               end do
         endif
      end do

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                       pf%biochem_pmin_col(c)
      end do


      ! column loop
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level phosphorus storage, for mass conservation check
!         col_endpb(c) = totpftp(c)
!         col_endpb(c) = sminp(c)
         col_endpb(c) = totcolp(c)
!         col_endpb(c) = totsomp(c)+cwdp(c)+totlitp(c)


         ! calculate total column-level inputs
         col_pinputs(c) = primp_to_labilep(c) + supplement_to_sminp(c) + supplement_to_sminp_surf_col(c)
!         col_pinputs(c) = leafp_to_litter_col(c)+frootp_to_litter_col(c)
!         col_pinputs(c) = sminp_to_plant(c)
!         col_pinputs(c) =  flux_mineralization_col(c)/dt
!         if (crop_prog) col_pinputs(c) = col_pinputs(c) + fert_p_to_sminp(c) 

         do p = col_pp%pfti(c), col_pp%pftf(c)
            if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
                col_pinputs(c) = col_pinputs(c) + supplement_to_plantp(p) * veg_pp%wtcol(p)
            end if
         end do

         ! forest fertilization
         if (forest_fert_exp) then
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 2 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               col_pinputs(c) = col_pinputs(c) + fert_dose(c,kmo)/dt
             end if
         end if

         if ((nu_com .ne. 'RD') .and. ECA_Pconst_RGspin) then
            col_poutputs(c) = secondp_to_occlp(c) + sminp_leached(c) + col_fire_ploss(c) + dwt_ploss(c) + product_ploss(c)
            do j = 1, nlevdecomp               
               col_poutputs(c) = col_poutputs(c) + &
                  (ps%solutionp_vr_col_cur(c,j) -  ps%solutionp_vr_col_prev(c,j)  + &
                  ps%labilep_vr_col_cur(c,j) -  ps%labilep_vr_col_prev(c,j) + &
                  ps%secondp_vr_col_cur(c,j) - ps%secondp_vr_col_prev(c,j) ) * dzsoi_decomp(j)/dt
            end do 
         else
            col_poutputs(c) = secondp_to_occlp(c) + sminp_leached(c) + col_fire_ploss(c) + dwt_ploss(c) + product_ploss(c)
         end if
         col_poutputs(c) = col_poutputs(c) +  sminp_runoff_col(c)
!         col_poutputs(c) = leafp_to_litter_col(c)+frootp_to_litter_col(c)
!         col_poutputs(c) =  flux_mineralization_col(c)/dt
!          col_poutputs(c) = sminp_to_plant(c)
         col_poutputs(c) = col_poutputs(c) - som_p_leached(c) + som_p_runoff_col(c)

!         if(is_active_betr_bgc .and. suplphos=='ALL')then
!           err_blg= phosphorusstate_vars%begblgp_col(c)- phosphorusstate_vars%totblgp_col(c) &
!             + primp_to_labilep(c) * dt &
!             + supplement_to_sminp(c) * dt &
!             + phosphorusflux_vars%pflx_plant_to_soilbgc_col(c)*dt &
!             - secondp_to_occlp(c) * dt  &
!             - sminp_leached(c) * dt + som_p_leached(c) * dt  &
!             - som_p_runoff_col(c)*dt-phosphorusflux_vars%sminp_to_plant_col(c)*dt &
!             - phosphorusflux_vars%fire_decomp_ploss_col(c)*dt
!           if(err_blg<0._r8)then
!             supplement_to_sminp(c) = supplement_to_sminp(c) - err_blg/dt
!             col_pinputs(c) = col_pinputs(c)  - err_blg/dt
!           endif
!         endif
         ! calculate the total column-level phosphorus balance error for this time step
         col_errpb(c) = (col_pinputs(c) - col_poutputs(c))*dt - &
              (col_endpb(c) - col_begpb(c))
!         if(c==15986  .and. get_nstep()== 9)then
!            err_found = .true.
!            err_index = c
!            exit
!         endif
         if (abs(col_errpb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
            exit
         end if
      end do ! end of columns loop


      if (err_found) then
         c = err_index
         write(iulog,*)'column pbalance error = ', col_errpb(c), c,  get_nstep()
         write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
         write(iulog,*)'begpb          = ',col_begpb(c)
         write(iulog,*)'endpb          = ',col_endpb(c)
         write(iulog,*)'begabgp        = ',phosphorusstate_vars%begabgp_col(c)
         write(iulog,*)'begblgp        = ',phosphorusstate_vars%begblgp_col(c)
         write(iulog,*)'delta store    = ',col_endpb(c)-col_begpb(c)
         write(iulog,*)'input mass     = ',col_pinputs(c)*dt
         write(iulog,*)'output mass    = ',col_poutputs(c)*dt
         write(iulog,*)'net flux       = ',(col_pinputs(c)-col_poutputs(c))*dt
         write(iulog,*)'weathering p   = ', primp_to_labilep(c) * dt
         write(iulog,*)'supp p         = ', supplement_to_sminp(c) * dt
         write(iulog,*)'supp p abv     = ', supplement_to_sminp_surf_col(c) * dt
         write(iulog,*)'2nd2occl p     = ', secondp_to_occlp(c) * dt
         write(iulog,*)'leached min p  = ',sminp_leached(c) * dt
         write(iulog,*)'lost to fire   = ', col_fire_ploss(c) * dt
         write(iulog,*)'harvest loss   = ', dwt_ploss(c) * dt
         write(iulog,*)'woodpro loss   = ', product_ploss(c) * dt
         write(iulog,*)'plant2bgcp     = ', phosphorusflux_vars%pflx_plant_to_soilbgc_col(c)*dt
         write(iulog,*)'leached som_p  = ', -som_p_leached(c) * dt
         write(iulog,*)'runoff  som_p  = ', som_p_runoff_col(c)*dt
         write(iulog,*)'soil2plant     = ', phosphorusflux_vars%sminp_to_plant_col(c)*dt
         write(iulog,*)'soil2plant_tran= ', phosphorusflux_vars%sminp_to_plant_trans_col(c)*dt
         write(iulog,*)'abgp           = ', phosphorusstate_vars%totabgp_col(c)
         write(iulog,*)'blgp           = ', phosphorusstate_vars%totblgp_col(c)
         write(iulog,*)'sminp_to_ppool = ', phosphorusflux_vars%sminp_to_ppool_col(c)*dt

         err_blg= phosphorusstate_vars%begblgp_col(c)- phosphorusstate_vars%totblgp_col(c) &
           + primp_to_labilep(c) * dt &
           + supplement_to_sminp(c) * dt &
           + phosphorusflux_vars%pflx_plant_to_soilbgc_col(c)*dt &
           - secondp_to_occlp(c) * dt  &
           - sminp_leached(c) * dt + som_p_leached(c) * dt  &
           - som_p_runoff_col(c)*dt-phosphorusflux_vars%sminp_to_plant_col(c)*dt &
           - phosphorusflux_vars%fire_decomp_ploss_col(c)*dt
         write(iulog,*)'err_blg =',err_blg

         err_abg=phosphorusstate_vars%begabgp_col(c)-phosphorusstate_vars%totabgp_col(c) &
           + phosphorusflux_vars%sminp_to_plant_col(c)*dt &
           + supplement_to_sminp_surf_col(c) * dt + &
           - phosphorusflux_vars%pflx_plant_to_soilbgc_col(c)*dt  &
           - dwt_ploss(c) * dt   &
           - product_ploss(c) * dt &
           - phosphorusflux_vars%fire_ploss_p2c_col(c) * dt
         write(iulog,*)'err_abg=',err_abg
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine PBalanceCheck

end module CNBalanceCheckMod
