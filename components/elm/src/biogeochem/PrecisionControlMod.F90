module PrecisionControlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use elm_varctl          , only : nu_com
  use elm_varpar          , only : ndecomp_pools
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
  subroutine PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, force leaf and deadstem c and n to 0 if
    ! they get too small.
    !
    ! !USES:
      !$acc routine seq
    use elm_varctl , only : iulog, use_c13, use_c14, use_fates
    use elm_varpar , only : nlevdecomp_full, crop_prog
    use pftvarcon  , only : crop, generic_crop, percrop
    use tracer_varcon          , only : is_active_betr_bgc
    use CNDecompCascadeConType , only : decomp_cascade_con
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
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
            if (abs(veg_cs%leafc(p)) < ccrit) then
               pc = pc + veg_cs%leafc(p)
               veg_cs%leafc(p) = 0._r8
               pn = pn + veg_ns%leafn(p)
               veg_ns%leafn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%leafc(p)
                  c13_veg_cs%leafc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%leafc(p)
                  c14_veg_cs%leafc(p) = 0._r8
               endif

               pp = pp + veg_ps%leafp(p)
               veg_ps%leafp(p) = 0._r8
            end if

            ! leaf storage C and N
            if (abs(veg_cs%leafc_storage(p)) < ccrit) then
               pc = pc + veg_cs%leafc_storage(p)
               veg_cs%leafc_storage(p) = 0._r8
               pn = pn + veg_ns%leafn_storage(p)
               veg_ns%leafn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%leafc_storage(p)
                  c13_veg_cs%leafc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%leafc_storage(p)
                  c14_veg_cs%leafc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%leafp_storage(p)
               veg_ps%leafp_storage(p) = 0._r8
            end if

            ! leaf transfer C and N
            if (abs(veg_cs%leafc_xfer(p)) < ccrit) then
               pc = pc + veg_cs%leafc_xfer(p)
               veg_cs%leafc_xfer(p) = 0._r8
               pn = pn + veg_ns%leafn_xfer(p)
               veg_ns%leafn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%leafc_xfer(p)
                  c13_veg_cs%leafc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%leafc_xfer(p)
                  c14_veg_cs%leafc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%leafp_xfer(p)
               veg_ps%leafp_xfer(p) = 0._r8
            end if

            ! froot C and N
            if (abs(veg_cs%frootc(p)) < ccrit) then
               pc = pc + veg_cs%frootc(p)
               veg_cs%frootc(p) = 0._r8
               pn = pn + veg_ns%frootn(p)
               veg_ns%frootn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%frootc(p)
                  c13_veg_cs%frootc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%frootc(p)
                  c14_veg_cs%frootc(p) = 0._r8
               endif

               pp = pp + veg_ps%frootp(p)
               veg_ps%frootp(p) = 0._r8
            end if

            ! froot storage C and N
            if (abs(veg_cs%frootc_storage(p)) < ccrit) then
               pc = pc + veg_cs%frootc_storage(p)
               veg_cs%frootc_storage(p) = 0._r8
               pn = pn + veg_ns%frootn_storage(p)
               veg_ns%frootn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%frootc_storage(p)
                  c13_veg_cs%frootc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%frootc_storage(p)
                  c14_veg_cs%frootc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%frootp_storage(p)
               veg_ps%frootp_storage(p) = 0._r8
            end if

            ! froot transfer C and N
            if (abs(veg_cs%frootc_xfer(p)) < ccrit) then
               pc = pc + veg_cs%frootc_xfer(p)
               veg_cs%frootc_xfer(p) = 0._r8
               pn = pn + veg_ns%frootn_xfer(p)
               veg_ns%frootn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%frootc_xfer(p)
                  c13_veg_cs%frootc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%frootc_xfer(p)
                  c14_veg_cs%frootc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%frootp_xfer(p)
               veg_ps%frootp_xfer(p) = 0._r8
            end if

            if ( crop_prog .and. &
                (crop(veg_pp%itype(p)) >= 1 .or. &
                 percrop(veg_pp%itype(p)) >= 1) )then
               ! grain C and N
               if (abs(veg_cs%grainc(p)) < ccrit) then
                  pc = pc + veg_cs%grainc(p)
                  veg_cs%grainc(p) = 0._r8
                  pn = pn + veg_ns%grainn(p)
                  veg_ns%grainn(p) = 0._r8
                  pp = pp + veg_ps%grainp(p)
                  veg_ps%grainp(p) = 0._r8
               end if

               ! grain storage C and N
               if (abs(veg_cs%grainc_storage(p)) < ccrit) then
                  pc = pc + veg_cs%grainc_storage(p)
                  veg_cs%grainc_storage(p) = 0._r8
                  pn = pn + veg_ns%grainn_storage(p)
                  veg_ns%grainn_storage(p) = 0._r8
                  pp = pp + veg_ps%grainp_storage(p)
                  veg_ps%grainp_storage(p) = 0._r8
               end if

               ! grain transfer C and N
               if (abs(veg_cs%grainc_xfer(p)) < ccrit) then
                  pc = pc + veg_cs%grainc_xfer(p)
                  veg_cs%grainc_xfer(p) = 0._r8
                  pn = pn + veg_ns%grainn_xfer(p)
                  veg_ns%grainn_xfer(p) = 0._r8
                  pp = pp + veg_ps%grainp_xfer(p)
                  veg_ps%grainp_xfer(p) = 0._r8
               end if
            end if

            ! livestem C and N
            if (abs(veg_cs%livestemc(p)) < ccrit) then
               pc = pc + veg_cs%livestemc(p)
               veg_cs%livestemc(p) = 0._r8
               pn = pn + veg_ns%livestemn(p)
               veg_ns%livestemn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%livestemc(p)
                  c13_veg_cs%livestemc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%livestemc(p)
                  c14_veg_cs%livestemc(p) = 0._r8
               endif

               pp = pp + veg_ps%livestemp(p)
               veg_ps%livestemp(p) = 0._r8
            end if

            ! livestem storage C and N
            if (abs(veg_cs%livestemc_storage(p)) < ccrit) then
               pc = pc + veg_cs%livestemc_storage(p)
               veg_cs%livestemc_storage(p) = 0._r8
               pn = pn + veg_ns%livestemn_storage(p)
               veg_ns%livestemn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%livestemc_storage(p)
                  c13_veg_cs%livestemc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%livestemc_storage(p)
                  c14_veg_cs%livestemc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%livestemp_storage(p)
               veg_ps%livestemp_storage(p) = 0._r8
            end if

            ! livestem transfer C and N
            if (abs(veg_cs%livestemc_xfer(p)) < ccrit) then
               pc = pc + veg_cs%livestemc_xfer(p)
               veg_cs%livestemc_xfer(p) = 0._r8
               pn = pn + veg_ns%livestemn_xfer(p)
               veg_ns%livestemn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%livestemc_xfer(p)
                  c13_veg_cs%livestemc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%livestemc_xfer(p)
                  c14_veg_cs%livestemc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%livestemp_xfer(p)
               veg_ps%livestemp_xfer(p) = 0._r8
            end if

            ! deadstem C and N
            if (abs(veg_cs%deadstemc(p)) < ccrit) then
               pc = pc + veg_cs%deadstemc(p)
               veg_cs%deadstemc(p) = 0._r8
               pn = pn + veg_ns%deadstemn(p)
               veg_ns%deadstemn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%deadstemc(p)
                  c13_veg_cs%deadstemc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%deadstemc(p)
                  c14_veg_cs%deadstemc(p) = 0._r8
               endif

               pp = pp + veg_ps%deadstemp(p)
               veg_ps%deadstemp(p) = 0._r8
            end if

            ! deadstem storage C and N
            if (abs(veg_cs%deadstemc_storage(p)) < ccrit) then
               pc = pc + veg_cs%deadstemc_storage(p)
               veg_cs%deadstemc_storage(p) = 0._r8
               pn = pn + veg_ns%deadstemn_storage(p)
               veg_ns%deadstemn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%deadstemc_storage(p)
                  c13_veg_cs%deadstemc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%deadstemc_storage(p)
                  c14_veg_cs%deadstemc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%deadstemp_storage(p)
               veg_ps%deadstemp_storage(p) = 0._r8
            end if

            ! deadstem transfer C and N
            if (abs(veg_cs%deadstemc_xfer(p)) < ccrit) then
               pc = pc + veg_cs%deadstemc_xfer(p)
               veg_cs%deadstemc_xfer(p) = 0._r8
               pn = pn + veg_ns%deadstemn_xfer(p)
               veg_ns%deadstemn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%deadstemc_xfer(p)
                  c13_veg_cs%deadstemc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%deadstemc_xfer(p)
                  c14_veg_cs%deadstemc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%deadstemp_xfer(p)
               veg_ps%deadstemp_xfer(p) = 0._r8
            end if

            ! livecroot C and N
            if (abs(veg_cs%livecrootc(p)) < ccrit) then
               pc = pc + veg_cs%livecrootc(p)
               veg_cs%livecrootc(p) = 0._r8
               pn = pn + veg_ns%livecrootn(p)
               veg_ns%livecrootn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%livecrootc(p)
                  c13_veg_cs%livecrootc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%livecrootc(p)
                  c14_veg_cs%livecrootc(p) = 0._r8
               endif

               pp = pp + veg_ps%livecrootp(p)
               veg_ps%livecrootp(p) = 0._r8
            end if

            ! livecroot storage C and N
            if (abs(veg_cs%livecrootc_storage(p)) < ccrit) then
               pc = pc + veg_cs%livecrootc_storage(p)
               veg_cs%livecrootc_storage(p) = 0._r8
               pn = pn + veg_ns%livecrootn_storage(p)
               veg_ns%livecrootn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%livecrootc_storage(p)
                  c13_veg_cs%livecrootc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%livecrootc_storage(p)
                  c14_veg_cs%livecrootc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%livecrootp_storage(p)
               veg_ps%livecrootp_storage(p) = 0._r8
            end if

            ! livecroot transfer C and N
            if (abs(veg_cs%livecrootc_xfer(p)) < ccrit) then
               pc = pc + veg_cs%livecrootc_xfer(p)
               veg_cs%livecrootc_xfer(p) = 0._r8
               pn = pn + veg_ns%livecrootn_xfer(p)
               veg_ns%livecrootn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%livecrootc_xfer(p)
                  c13_veg_cs%livecrootc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%livecrootc_xfer(p)
                  c14_veg_cs%livecrootc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%livecrootp_xfer(p)
               veg_ps%livecrootp_xfer(p) = 0._r8
            end if

            ! deadcroot C and N
            if (abs(veg_cs%deadcrootc(p)) < ccrit) then
               pc = pc + veg_cs%deadcrootc(p)
               veg_cs%deadcrootc(p) = 0._r8
               pn = pn + veg_ns%deadcrootn(p)
               veg_ns%deadcrootn(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%deadcrootc(p)
                  c13_veg_cs%deadcrootc(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%deadcrootc(p)
                  c14_veg_cs%deadcrootc(p) = 0._r8
               endif

               pp = pp + veg_ps%deadcrootp(p)
               veg_ps%deadcrootp(p) = 0._r8
            end if

            ! deadcroot storage C and N
            if (abs(veg_cs%deadcrootc_storage(p)) < ccrit) then
               pc = pc + veg_cs%deadcrootc_storage(p)
               veg_cs%deadcrootc_storage(p) = 0._r8
               pn = pn + veg_ns%deadcrootn_storage(p)
               veg_ns%deadcrootn_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%deadcrootc_storage(p)
                  c13_veg_cs%deadcrootc_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%deadcrootc_storage(p)
                  c14_veg_cs%deadcrootc_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%deadcrootp_storage(p)
               veg_ps%deadcrootp_storage(p) = 0._r8
            end if

            ! deadcroot transfer C and N
            if (abs(veg_cs%deadcrootc_xfer(p)) < ccrit) then
               pc = pc + veg_cs%deadcrootc_xfer(p)
               veg_cs%deadcrootc_xfer(p) = 0._r8
               pn = pn + veg_ns%deadcrootn_xfer(p)
               veg_ns%deadcrootn_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%deadcrootc_xfer(p)
                  c13_veg_cs%deadcrootc_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%deadcrootc_xfer(p)
                  c14_veg_cs%deadcrootc_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%deadcrootp_xfer(p)
               veg_ps%deadcrootp_xfer(p) = 0._r8
            end if

            ! gresp_storage (C only)
            if (abs(veg_cs%gresp_storage(p)) < ccrit) then
               pc = pc + veg_cs%gresp_storage(p)
               veg_cs%gresp_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%gresp_storage(p)
                  c13_veg_cs%gresp_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%gresp_storage(p)
                  c14_veg_cs%gresp_storage(p) = 0._r8
               endif
            end if

            ! gresp_xfer(c only)
            if (abs(veg_cs%gresp_xfer(p)) < ccrit) then
               pc = pc + veg_cs%gresp_xfer(p)
               veg_cs%gresp_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%gresp_xfer(p)
                  c13_veg_cs%gresp_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%gresp_xfer(p)
                  c14_veg_cs%gresp_xfer(p) = 0._r8
               endif
            end if

            ! cpool (C only)
            if (abs(veg_cs%cpool(p)) < ccrit) then
               pc = pc + veg_cs%cpool(p)
               veg_cs%cpool(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%cpool(p)
                  c13_veg_cs%cpool(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%cpool(p)
                  c14_veg_cs%cpool(p) = 0._r8
               endif
            end if

            if ( crop_prog .and. &
                (crop(veg_pp%itype(p)) >= 1 .or. &
                 percrop(veg_pp%itype(p)) >= 1) )then
               ! xsmrpool (C only)
               if (abs(veg_cs%xsmrpool(p)) < ccrit) then
                  pc = pc + veg_cs%xsmrpool(p)
                  veg_cs%xsmrpool(p) = 0._r8
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

            veg_cs%ctrunc(p) = veg_cs%ctrunc(p) + pc
            veg_ns%ntrunc(p) = veg_ns%ntrunc(p) + pn
            veg_ps%ptrunc(p) = veg_ps%ptrunc(p) + pp

            if ( use_c13 ) then
               c13_veg_cs%ctrunc(p) = c13_veg_cs%ctrunc(p) + pc13
            endif
            if ( use_c14 ) then
               c14_veg_cs%ctrunc(p) = c14_veg_cs%ctrunc(p) + pc14
            endif

         end do ! end of pft loop
      end if ! end of if(.not.use_fates)

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

                  if (abs(col_cs%decomp_cpools_vr(c,j,k)) < ccrit) then
                     cc = cc + col_cs%decomp_cpools_vr(c,j,k)
                     col_cs%decomp_cpools_vr(c,j,k) = 0._r8
                     cn = cn + col_ns%decomp_npools_vr(c,j,k)
                     col_ns%decomp_npools_vr(c,j,k) = 0._r8
                     if ( use_c13 ) then
                        cc13 = cc13 + c13_col_cs%decomp_cpools_vr(c,j,k)
                        c13_col_cs%decomp_cpools_vr(c,j,k) = 0._r8
                     endif
                     if ( use_c14 ) then
                        cc14 = cc14 + c14_col_cs%decomp_cpools_vr(c,j,k)
                        c14_col_cs%decomp_cpools_vr(c,j,k) = 0._r8
                     endif
                  end if

               end do

               ! not doing precision control on soil mineral N, since it will
               ! be getting the N truncation flux anyway.

               col_cs%ctrunc_vr(c,j) = col_cs%ctrunc_vr(c,j) + cc
               col_ns%ntrunc_vr(c,j) = col_ns%ntrunc_vr(c,j) + cn
               if ( use_c13 ) then
                  c13_col_cs%ctrunc_vr(c,j) = c13_col_cs%ctrunc_vr(c,j) + cc13
               endif
               if ( use_c14 ) then
                  c14_col_cs%ctrunc_vr(c,j) = c14_col_cs%ctrunc_vr(c,j) + cc14
               endif
            end do

         end do   ! end of column loop

         ! remove small negative perturbations for stability purposes, if any should arise.
         
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            do j = 1,nlevdecomp_full
               if (abs(col_ns%smin_no3_vr(c,j)) < ncrit/1e4_r8) then
                  if ( col_ns%smin_no3_vr(c,j)  < 0._r8 ) then
#ifndef _OPENACC
                     write(iulog, *) '-10^-12 < smin_no3 < 0. resetting to zero.'
                     write(iulog, *) 'smin_no3_vr_col(c,j), c, j: ', col_ns%smin_no3_vr(c,j), c, j
                     col_ns%smin_no3_vr(c,j) = 0._r8
#endif
                  endif
               end if
               if (abs(col_ns%smin_nh4_vr(c,j)) < ncrit/1e4_r8) then
                  if ( col_ns%smin_nh4_vr(c,j)  < 0._r8 ) then
#ifndef _OPENACC
                     write(iulog, *) '-10^-12 < smin_nh4 < 0. resetting to zero.'
                     write(iulog, *) 'smin_nh4_vr_col(c,j), c, j: ', col_ns%smin_nh4_vr(c,j), c, j
                     col_ns%smin_nh4_vr(c,j) = 0._r8
#endif
                  endif
               end if
            end do
         end do


        if (nu_com .eq. 'ECA') then
            ! decompose P pool adjust according to C pool
            !do fc = 1,num_soilc
            !   c = filter_soilc(fc)
            !   do j = 1,nlevdecomp_full
            !      cp_eca = 0.0_r8
            !      do l = 1,ndecomp_pools
            !         if (abs(col_cs%decomp_cpools_vr(c,j,k)) < ccrit) then
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
                     if ( col_cs%decomp_cpools_vr(c,j,l) > 0.0_r8 ) then
                          if(abs(col_cs%decomp_cpools_vr(c,j,l) / col_ns%decomp_npools_vr(c,j,l) - initial_cn_ratio(l) ) > 1.0e-3_r8 &
                          .and. (.not. floating_cn_ratio_decomp_pools(l)) ) then
                        cn_eca = cn_eca - ( col_cs%decomp_cpools_vr(c,j,l) / initial_cn_ratio(l) - col_ns%decomp_npools_vr(c,j,l) )

                          col_ns%decomp_npools_vr(c,j,l) = col_cs%decomp_cpools_vr(c,j,l) / initial_cn_ratio(l)
                     end if
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
#ifndef _OPENACC                                
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_npools is negative: ',j,l,col_ns%decomp_npools_vr(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
                        end if
                     end if
                     if ( col_ps%decomp_ppools_vr(c,j,l)  < 0.0_r8 .and. floating_cp_ratio_decomp_pools(l) ) then
                        if ( abs(col_ps%decomp_ppools_vr(c,j,l))  < ncrit/1e4_r8 ) then
                           cp_eca = cp_eca - ncrit/1e4_r8 + col_ps%decomp_ppools_vr(c,j,l)
                           col_ps%decomp_ppools_vr(c,j,l) = ncrit/1e4_r8
                         else
#ifndef _OPENACC
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_ppools is negative: ',j,l,col_ps%decomp_ppools_vr(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
                         end if
                     end if

                  end do

                  col_ns%ntrunc_vr(c,j) = col_ns%ntrunc_vr(c,j) + cn_eca
                  col_ps%ptrunc_vr(c,j) = col_ps%ptrunc_vr(c,j) + cp_eca

               end do
            end do

          if(.not.use_fates) then
            do fp = 1,num_soilp
               p = filter_soilp(fp)
#ifndef _OPENACC
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
#endif
            end do
         endif

       endif  !if ECA

      endif ! if (.not. is_active_betr_bgc)


    end associate

 end subroutine PrecisionControl

end module PrecisionControlMod
