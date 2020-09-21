module PrecisionControlMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables 
  ! 
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use clm_varctl          , only : nu_com
  use elm_varpar          , only : ndecomp_pools
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PhosphorusStateType , only : phosphorusstate_type
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType      , only : col_ns, col_ps
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType  , only : veg_ns, veg_ps

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: PrecisionControl
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars,&
       phosphorusstate_vars)
    !
    ! !DESCRIPTION: 
    ! On the radiation time step, force leaf and deadstem c and n to 0 if
    ! they get too small.
    !
    ! !USES:
    use clm_varctl , only : iulog, use_c13, use_c14, use_nitrif_denitrif, use_fates
    use elm_varpar , only : nlevdecomp_full, crop_prog
    use pftvarcon  , only : nc3crop
    use tracer_varcon          , only : is_active_betr_bgc    
    use CNDecompCascadeConType , only : decomp_cascade_con
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l  ! indices
    integer :: fp,fc    ! lake filter indices
    real(r8):: pc,pn,pp    ! truncation terms for patch-level corrections
    real(r8):: cc,cn,cp    ! truncation terms for column-level corrections
    real(r8):: pc13     ! truncation terms for patch-level corrections
    real(r8):: cc13     ! truncation terms for column-level corrections
    real(r8):: pc14     ! truncation terms for patch-level corrections
    real(r8):: cc14     ! truncation terms for column-level corrections
    real(r8):: ccrit    ! critical carbon state value for truncation
    real(r8):: ncrit    ! critical nitrogen state value for truncation
    real(r8):: pcrit    ! critical phosphorus state value for truncation
    real(r8):: cc_eca
    real(r8):: cn_eca
    real(r8):: cp_eca
    !-----------------------------------------------------------------------

    associate(&
         cs      => carbonstate_vars     , &
         csv2    => col_cs               , &
         vcsv2   => veg_cs               , &
         ns      => nitrogenstate_vars   , &
         ps      => phosphorusstate_vars , &
         c13cs   => c13_carbonstate_vars , &
         c13csv2 => c13_col_cs           , &
         c13vcsv2=> c13_veg_cs           , &
         c14cs   => c14_carbonstate_vars , &
         c14csv2 => c14_col_cs           , &
         c14vcsv2=> c14_veg_cs           , &
         floating_cn_ratio_decomp_pools   =>    decomp_cascade_con%floating_cn_ratio_decomp_pools , &
         floating_cp_ratio_decomp_pools   =>    decomp_cascade_con%floating_cp_ratio_decomp_pools , &
         initial_cn_ratio                 =>    decomp_cascade_con%initial_cn_ratio                 &
         )

      ! set the critical carbon state value for truncation (gC/m2)
      ccrit = 1.e-8_r8

      ! set the critical nitrogen state value for truncation (gN/m2)
      ncrit = 1.e-8_r8

      ! set the critical phosphorus state value for truncation (gN/m2)
      pcrit = 1.e-8_r8

      ! patch loop
      if (.not.use_fates) then
         do fp = 1,num_soilp
            p = filter_soilp(fp)

            ! initialize the patch-level C and N truncation terms
            pc = 0._r8
            pn = 0._r8
            pp = 0._r8
            if ( use_c13 ) pc13 = 0._r8
            if ( use_c14 ) pc14 = 0._r8

            ! do tests on state variables for precision control
            ! for linked C-N state variables, perform precision test on
            ! the C component, but truncate C, C13, and N components

            ! leaf C and N
            if (abs(vcsv2%leafc(p)) < ccrit) then
               pc = pc + vcsv2%leafc(p)
               vcsv2%leafc(p) = 0._r8
               pn = pn + veg_ns%leafn(p)
               veg_ns%leafn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%leafc(p)
                  c13vcsv2%leafc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%leafc(p)
                  c14vcsv2%leafc(p) = 0._r8
               endif

               pp = pp + veg_ps%leafp(p)
               veg_ps%leafp(p) = 0._r8
            end if

            ! leaf storage C and N
            if (abs(vcsv2%leafc_storage(p)) < ccrit) then
               pc = pc + vcsv2%leafc_storage(p)
               vcsv2%leafc_storage(p) = 0._r8
               pn = pn + veg_ns%leafn_storage(p)
               veg_ns%leafn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%leafc_storage(p)
                  c13vcsv2%leafc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%leafc_storage(p)
                  c14vcsv2%leafc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%leafp_storage(p)
               veg_ps%leafp_storage(p) = 0._r8
            end if

            ! leaf transfer C and N
            if (abs(vcsv2%leafc_xfer(p)) < ccrit) then
               pc = pc + vcsv2%leafc_xfer(p)
               vcsv2%leafc_xfer(p) = 0._r8
               pn = pn + veg_ns%leafn_xfer(p)
               veg_ns%leafn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%leafc_xfer(p)
                  c13vcsv2%leafc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%leafc_xfer(p)
                  c14vcsv2%leafc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%leafp_xfer(p)
               veg_ps%leafp_xfer(p) = 0._r8
            end if

            ! froot C and N
            if (abs(vcsv2%frootc(p)) < ccrit) then
               pc = pc + vcsv2%frootc(p)
               vcsv2%frootc(p) = 0._r8
               pn = pn + veg_ns%frootn(p)
               veg_ns%frootn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%frootc(p)
                  c13vcsv2%frootc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%frootc(p)
                  c14vcsv2%frootc(p) = 0._r8
               endif

               pp = pp + veg_ps%frootp(p)
               veg_ps%frootp(p) = 0._r8
            end if

            ! froot storage C and N
            if (abs(vcsv2%frootc_storage(p)) < ccrit) then
               pc = pc + vcsv2%frootc_storage(p)
               vcsv2%frootc_storage(p) = 0._r8
               pn = pn + veg_ns%frootn_storage(p)
               veg_ns%frootn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%frootc_storage(p)
                  c13vcsv2%frootc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%frootc_storage(p)
                  c14vcsv2%frootc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%frootp_storage(p)
               veg_ps%frootp_storage(p) = 0._r8
            end if

            ! froot transfer C and N
            if (abs(vcsv2%frootc_xfer(p)) < ccrit) then
               pc = pc + vcsv2%frootc_xfer(p)
               vcsv2%frootc_xfer(p) = 0._r8
               pn = pn + veg_ns%frootn_xfer(p)
               veg_ns%frootn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%frootc_xfer(p)
                  c13vcsv2%frootc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%frootc_xfer(p)
                  c14vcsv2%frootc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%frootp_xfer(p)
               veg_ps%frootp_xfer(p) = 0._r8
            end if

            if ( crop_prog .and. veg_pp%itype(p) >= nc3crop )then
               ! grain C and N
               if (abs(vcsv2%grainc(p)) < ccrit) then
                  pc = pc + vcsv2%grainc(p)
                  vcsv2%grainc(p) = 0._r8
                  pn = pn + veg_ns%grainn(p)
                  veg_ns%grainn(p) = 0._r8
                  pp = pp + veg_ps%grainp(p)
                  veg_ps%grainp(p) = 0._r8
               end if

               ! grain storage C and N
               if (abs(vcsv2%grainc_storage(p)) < ccrit) then
                  pc = pc + vcsv2%grainc_storage(p)
                  vcsv2%grainc_storage(p) = 0._r8
                  pn = pn + veg_ns%grainn_storage(p)
                  veg_ns%grainn_storage(p) = 0._r8
                  pp = pp + veg_ps%grainp_storage(p)
                  veg_ps%grainp_storage(p) = 0._r8
               end if

               ! grain transfer C and N
               if (abs(vcsv2%grainc_xfer(p)) < ccrit) then
                  pc = pc + vcsv2%grainc_xfer(p)
                  vcsv2%grainc_xfer(p) = 0._r8
                  pn = pn + veg_ns%grainn_xfer(p)
                  veg_ns%grainn_xfer(p) = 0._r8
                  pp = pp + veg_ps%grainp_xfer(p)
                  veg_ps%grainp_xfer(p) = 0._r8
               end if
            end if

            ! livestem C and N
            if (abs(vcsv2%livestemc(p)) < ccrit) then
               pc = pc + vcsv2%livestemc(p)
               vcsv2%livestemc(p) = 0._r8
               pn = pn + veg_ns%livestemn(p)
               veg_ns%livestemn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%livestemc(p)
                  c13vcsv2%livestemc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%livestemc(p)
                  c14vcsv2%livestemc(p) = 0._r8
               endif

               pp = pp + veg_ps%livestemp(p)
               veg_ps%livestemp(p) = 0._r8
            end if

            ! livestem storage C and N
            if (abs(vcsv2%livestemc_storage(p)) < ccrit) then
               pc = pc + vcsv2%livestemc_storage(p)
               vcsv2%livestemc_storage(p) = 0._r8
               pn = pn + veg_ns%livestemn_storage(p)
               veg_ns%livestemn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%livestemc_storage(p)
                  c13vcsv2%livestemc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%livestemc_storage(p)
                  c14vcsv2%livestemc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%livestemp_storage(p)
               veg_ps%livestemp_storage(p) = 0._r8
            end if

            ! livestem transfer C and N
            if (abs(vcsv2%livestemc_xfer(p)) < ccrit) then
               pc = pc + vcsv2%livestemc_xfer(p)
               vcsv2%livestemc_xfer(p) = 0._r8
               pn = pn + veg_ns%livestemn_xfer(p)
               veg_ns%livestemn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%livestemc_xfer(p)
                  c13vcsv2%livestemc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%livestemc_xfer(p)
                  c14vcsv2%livestemc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%livestemp_xfer(p)
               veg_ps%livestemp_xfer(p) = 0._r8
            end if

            ! deadstem C and N
            if (abs(vcsv2%deadstemc(p)) < ccrit) then
               pc = pc + vcsv2%deadstemc(p)
               vcsv2%deadstemc(p) = 0._r8
               pn = pn + veg_ns%deadstemn(p)
               veg_ns%deadstemn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%deadstemc(p)
                  c13vcsv2%deadstemc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%deadstemc(p)
                  c14vcsv2%deadstemc(p) = 0._r8
               endif

               pp = pp + veg_ps%deadstemp(p)
               veg_ps%deadstemp(p) = 0._r8
            end if

            ! deadstem storage C and N
            if (abs(vcsv2%deadstemc_storage(p)) < ccrit) then
               pc = pc + vcsv2%deadstemc_storage(p)
               vcsv2%deadstemc_storage(p) = 0._r8
               pn = pn + veg_ns%deadstemn_storage(p)
               veg_ns%deadstemn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%deadstemc_storage(p)
                  c13vcsv2%deadstemc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%deadstemc_storage(p)
                  c14vcsv2%deadstemc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%deadstemp_storage(p)
               veg_ps%deadstemp_storage(p) = 0._r8
            end if

            ! deadstem transfer C and N
            if (abs(vcsv2%deadstemc_xfer(p)) < ccrit) then
               pc = pc + vcsv2%deadstemc_xfer(p)
               vcsv2%deadstemc_xfer(p) = 0._r8
               pn = pn + veg_ns%deadstemn_xfer(p)
               veg_ns%deadstemn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%deadstemc_xfer(p)
                  c13vcsv2%deadstemc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%deadstemc_xfer(p)
                  c14vcsv2%deadstemc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%deadstemp_xfer(p)
               veg_ps%deadstemp_xfer(p) = 0._r8
            end if

            ! livecroot C and N
            if (abs(vcsv2%livecrootc(p)) < ccrit) then
               pc = pc + vcsv2%livecrootc(p)
               vcsv2%livecrootc(p) = 0._r8
               pn = pn + veg_ns%livecrootn(p)
               veg_ns%livecrootn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%livecrootc(p)
                  c13vcsv2%livecrootc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%livecrootc(p)
                  c14vcsv2%livecrootc(p) = 0._r8
               endif

               pp = pp + veg_ps%livecrootp(p)
               veg_ps%livecrootp(p) = 0._r8
            end if

            ! livecroot storage C and N
            if (abs(vcsv2%livecrootc_storage(p)) < ccrit) then
               pc = pc + vcsv2%livecrootc_storage(p)
               vcsv2%livecrootc_storage(p) = 0._r8
               pn = pn + veg_ns%livecrootn_storage(p)
               veg_ns%livecrootn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%livecrootc_storage(p)
                  c13vcsv2%livecrootc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%livecrootc_storage(p)
                  c14vcsv2%livecrootc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%livecrootp_storage(p)
               veg_ps%livecrootp_storage(p) = 0._r8
            end if

            ! livecroot transfer C and N
            if (abs(vcsv2%livecrootc_xfer(p)) < ccrit) then
               pc = pc + vcsv2%livecrootc_xfer(p)
               vcsv2%livecrootc_xfer(p) = 0._r8
               pn = pn + veg_ns%livecrootn_xfer(p)
               veg_ns%livecrootn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%livecrootc_xfer(p)
                  c13vcsv2%livecrootc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%livecrootc_xfer(p)
                  c14vcsv2%livecrootc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%livecrootp_xfer(p)
               veg_ps%livecrootp_xfer(p) = 0._r8
            end if

            ! deadcroot C and N
            if (abs(vcsv2%deadcrootc(p)) < ccrit) then
               pc = pc + vcsv2%deadcrootc(p)
               vcsv2%deadcrootc(p) = 0._r8
               pn = pn + veg_ns%deadcrootn(p)
               veg_ns%deadcrootn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%deadcrootc(p)
                  c13vcsv2%deadcrootc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%deadcrootc(p)
                  c14vcsv2%deadcrootc(p) = 0._r8
               endif

               pp = pp + veg_ps%deadcrootp(p)
               veg_ps%deadcrootp(p) = 0._r8
            end if

            ! deadcroot storage C and N
            if (abs(vcsv2%deadcrootc_storage(p)) < ccrit) then
               pc = pc + vcsv2%deadcrootc_storage(p)
               vcsv2%deadcrootc_storage(p) = 0._r8
               pn = pn + veg_ns%deadcrootn_storage(p)
               veg_ns%deadcrootn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%deadcrootc_storage(p)
                  c13vcsv2%deadcrootc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%deadcrootc_storage(p)
                  c14vcsv2%deadcrootc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%deadcrootp_storage(p)
               veg_ps%deadcrootp_storage(p) = 0._r8
            end if

            ! deadcroot transfer C and N
            if (abs(vcsv2%deadcrootc_xfer(p)) < ccrit) then
               pc = pc + vcsv2%deadcrootc_xfer(p)
               vcsv2%deadcrootc_xfer(p) = 0._r8
               pn = pn + veg_ns%deadcrootn_xfer(p)
               veg_ns%deadcrootn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%deadcrootc_xfer(p)
                  c13vcsv2%deadcrootc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%deadcrootc_xfer(p)
                  c14vcsv2%deadcrootc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%deadcrootp_xfer(p)
               veg_ps%deadcrootp_xfer(p) = 0._r8
            end if

            ! gresp_storage (C only)
            if (abs(vcsv2%gresp_storage(p)) < ccrit) then
               pc = pc + vcsv2%gresp_storage(p)
               vcsv2%gresp_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%gresp_storage(p)
                  c13vcsv2%gresp_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%gresp_storage(p)
                  c14vcsv2%gresp_storage(p) = 0._r8
               endif
            end if

            ! gresp_xfer(c only)
            if (abs(vcsv2%gresp_xfer(p)) < ccrit) then
               pc = pc + vcsv2%gresp_xfer(p)
               vcsv2%gresp_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%gresp_xfer(p)
                  c13vcsv2%gresp_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%gresp_xfer(p)
                  c14vcsv2%gresp_xfer(p) = 0._r8
               endif
            end if

            ! cpool (C only)
            if (abs(vcsv2%cpool(p)) < ccrit) then
               pc = pc + vcsv2%cpool(p)
               vcsv2%cpool(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13vcsv2%cpool(p)
                  c13vcsv2%cpool(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14vcsv2%cpool(p)
                  c14vcsv2%cpool(p) = 0._r8
               endif
            end if

            if ( crop_prog .and. veg_pp%itype(p) >= nc3crop )then
               ! xsmrpool (C only)
               if (abs(vcsv2%xsmrpool(p)) < ccrit) then
                  pc = pc + vcsv2%xsmrpool(p)
                  vcsv2%xsmrpool(p) = 0._r8
               end if
            end if

            ! retransn (N only)
            if (abs(veg_ns%retransn(p)) < ncrit) then
               pn = pn + veg_ns%retransn(p)
               veg_ns%retransn(p) = 0._r8
            end if

            ! retransp (P only)
            if (abs(veg_ps%retransp(p)) < pcrit) then
               pp = pp + veg_ps%retransp(p)
               veg_ps%retransp(p) = 0._r8
            end if

            ! npool (N only)
            if (abs(veg_ns%npool(p)) < ncrit) then
               pn = pn + veg_ns%npool(p)
               veg_ns%npool(p) = 0._r8
            end if

            ! ppool (P only)
            if (abs(veg_ps%ppool(p)) < pcrit) then
               pp = pp + veg_ps%ppool(p)
               veg_ps%ppool(p) = 0._r8
            end if

            vcsv2%ctrunc(p) = vcsv2%ctrunc(p) + pc
            veg_ns%ntrunc(p) = veg_ns%ntrunc(p) + pn
            veg_ps%ptrunc(p) = veg_ps%ptrunc(p) + pp

            if ( use_c13 ) then
               c13vcsv2%ctrunc(p) = c13vcsv2%ctrunc(p) + pc13
            endif
            if ( use_c14 ) then
               c14vcsv2%ctrunc(p) = c14vcsv2%ctrunc(p) + pc14
            endif

         end do ! end of pft loop
      end if ! end of if(not.use_fates)

      if (.not. is_active_betr_bgc) then

         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            do j = 1,nlevdecomp_full
               ! initialize the column-level C and N truncation terms
               cc = 0._r8
               if ( use_c13 ) cc13 = 0._r8
               if ( use_c14 ) cc14 = 0._r8
               cn = 0._r8

               ! do tests on state variables for precision control
               ! for linked C-N state variables, perform precision test on
               ! the C component, but truncate both C and N components


               ! all decomposing pools C and N
               do k = 1, ndecomp_pools

                  if (abs(csv2%decomp_cpools_vr(c,j,k)) < ccrit) then
                     cc = cc + csv2%decomp_cpools_vr(c,j,k)
                     csv2%decomp_cpools_vr(c,j,k) = 0._r8
                     if (.not.use_fates) then
                        cn = cn + col_ns%decomp_npools_vr(c,j,k)
                        col_ns%decomp_npools_vr(c,j,k) = 0._r8
                     endif
                     if ( use_c13 ) then
                        cc13 = cc13 + c13csv2%decomp_cpools_vr(c,j,k)
                        c13csv2%decomp_cpools_vr(c,j,k) = 0._r8
                     endif
                     if ( use_c14 ) then
                        cc14 = cc14 + c14csv2%decomp_cpools_vr(c,j,k)
                        c14csv2%decomp_cpools_vr(c,j,k) = 0._r8
                     endif
                  end if

               end do

               ! not doing precision control on soil mineral N, since it will
               ! be getting the N truncation flux anyway.

               csv2%ctrunc_vr(c,j) = csv2%ctrunc_vr(c,j) + cc
               if (.not.use_fates) then
                  col_ns%ntrunc_vr(c,j) = col_ns%ntrunc_vr(c,j) + cn
               endif
               if ( use_c13 ) then
                  c13csv2%ctrunc_vr(c,j) = c13csv2%ctrunc_vr(c,j) + cc13
               endif
               if ( use_c14 ) then
                  c14csv2%ctrunc_vr(c,j) = c14csv2%ctrunc_vr(c,j) + cc14
               endif
            end do

         end do   ! end of column loop

         if (use_nitrif_denitrif) then
            ! remove small negative perturbations for stability purposes, if any should arise.

            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full
                  if (abs(col_ns%smin_no3_vr(c,j)) < ncrit/1e4_r8) then
                     if ( col_ns%smin_no3_vr(c,j)  < 0._r8 ) then
                        write(iulog, *) '-10^-12 < smin_no3 < 0. resetting to zero.'
                        write(iulog, *) 'smin_no3_vr_col(c,j), c, j: ', col_ns%smin_no3_vr(c,j), c, j
                        col_ns%smin_no3_vr(c,j) = 0._r8
                     endif
                  end if
                  if (abs(col_ns%smin_nh4_vr(c,j)) < ncrit/1e4_r8) then
                     if ( col_ns%smin_nh4_vr(c,j)  < 0._r8 ) then
                        write(iulog, *) '-10^-12 < smin_nh4 < 0. resetting to zero.'
                        write(iulog, *) 'smin_nh4_vr_col(c,j), c, j: ', col_ns%smin_nh4_vr(c,j), c, j
                        col_ns%smin_nh4_vr(c,j) = 0._r8
                     endif
                  end if
               end do
            end do
         endif

         if (nu_com .eq. 'ECA') then
            ! decompose P pool adjust according to C pool
            !do fc = 1,num_soilc
            !   c = filter_soilc(fc)
            !   do j = 1,nlevdecomp_full
            !      cp_eca = 0.0_r8
            !      do l = 1,ndecomp_pools
            !         if (abs(csv2%decomp_cpools_vr(c,j,k)) < ccrit) then
            !            if (.not.use_fates) then
            !               cp_eca = cp_eca + col_ps%decomp_ppools_vr(c,j,k)
            !               col_ps%decomp_ppools_vr(c,j,k) = 0._r8
            !            endif
            !         endif
            !      end do
            !      col_ps%ptrunc_vr(c,j) = col_ps%ptrunc_vr(c,j) + cp_eca
            !   end do
            !end do

            ! fix soil CN ratio drift (normally < 0.01% drift)
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full
                  cn_eca = 0.0_r8
                  do l = 1,ndecomp_pools
                     if ( csv2%decomp_cpools_vr(c,j,l) > 0.0_r8 .and.  &
                          abs(csv2%decomp_cpools_vr(c,j,l) / col_ns%decomp_npools_vr(c,j,l) - initial_cn_ratio(l) ) > 1.0e-3_r8 &
                          .and. (.not. floating_cn_ratio_decomp_pools(l)) ) then
                        cn_eca = cn_eca - ( csv2%decomp_cpools_vr(c,j,l) / initial_cn_ratio(l) - col_ns%decomp_npools_vr(c,j,l) )
                        col_ns%decomp_npools_vr(c,j,l) = csv2%decomp_cpools_vr(c,j,l) / initial_cn_ratio(l)
                     end if
                  end do
                  col_ns%ntrunc_vr(c,j) = col_ns%ntrunc_vr(c,j) + cn_eca
               end do
             end do

            ! remove small negative perturbations for stability purposes, if any should arise in N,P pools
            ! for floating CN, CP ratio pools
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full

                  cn_eca = 0.0_r8
                  cp_eca = 0.0_r8
                  do l = 1,ndecomp_pools
                     if ( col_ns%decomp_npools_vr(c,j,l) < 0.0_r8 .and. floating_cn_ratio_decomp_pools(l) ) then
                        if ( abs(col_ns%decomp_npools_vr(c,j,l))  < ncrit ) then
                           cn_eca = cn_eca - ncrit + col_ns%decomp_npools_vr(c,j,l)
                           col_ns%decomp_npools_vr(c,j,l) = ncrit
                        else
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_npools is negative: ',j,l,col_ns%decomp_npools_vr(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
                        end if
                     end if
                     if ( col_ps%decomp_ppools_vr(c,j,l)  < 0.0_r8 .and. floating_cp_ratio_decomp_pools(l) ) then
                        if ( abs(col_ps%decomp_ppools_vr(c,j,l))  < ncrit/1e4_r8 ) then
                           cp_eca = cp_eca - ncrit/1e4_r8 + col_ps%decomp_ppools_vr(c,j,l)
                           col_ps%decomp_ppools_vr(c,j,l) = ncrit/1e4_r8
                         else 
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_ppools is negative: ',j,l,col_ps%decomp_ppools_vr(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
                         end if
                     end if

                  end do

                  col_ns%ntrunc_vr(c,j) = col_ns%ntrunc_vr(c,j) + cn_eca
                  col_ps%ptrunc_vr(c,j) = col_ps%ptrunc_vr(c,j) + cp_eca

               end do
            end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if (veg_ns%retransn(p) < 0._r8) then
                  write(iulog, *) 'error retransn_patch is negative: ',p
                  write(iulog, *) 'retransn_patch: ', veg_ns%retransn(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (veg_ns%npool(p) < 0._r8) then
                  write(iulog, *) 'error npool_patch is negative: ',p
                  write(iulog, *) 'npool_patch: ', veg_ns%npool(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (veg_ps%retransp(p) < 0._r8) then
                  write(iulog, *) 'error retransp_patch is negative: ',p
                  write(iulog, *) 'retransp_patch: ', veg_ps%retransp(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (veg_ps%ppool(p) < 0._r8) then
                  write(iulog, *) 'error ppool_patch is negative: ',p
                  write(iulog, *) 'ppool_patch: ', veg_ps%ppool(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
            end do

         endif

      endif ! if (.not. is_active_betr_bgc)

    end associate

 end subroutine PrecisionControl

end module PrecisionControlMod
