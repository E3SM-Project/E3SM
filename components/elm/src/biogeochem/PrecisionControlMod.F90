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
    use pftvarcon  , only : iscft
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
            if (abs(veg_cs%core%leaf(p)) < ccrit) then
               pc = pc + veg_cs%core%leaf(p)
               veg_cs%core%leaf(p) = 0._r8
               pn = pn + veg_ns%core%leaf(p)
               veg_ns%core%leaf(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%leaf(p)
                  c13_veg_cs%core%leaf(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%leaf(p)
                  c14_veg_cs%core%leaf(p) = 0._r8
               endif

               pp = pp + veg_ps%core%leaf(p)
               veg_ps%core%leaf(p) = 0._r8
            end if

            ! leaf storage C and N
            if (abs(veg_cs%core%leaf_storage(p)) < ccrit) then
               pc = pc + veg_cs%core%leaf_storage(p)
               veg_cs%core%leaf_storage(p) = 0._r8
               pn = pn + veg_ns%core%leaf_storage(p)
               veg_ns%core%leaf_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%leaf_storage(p)
                  c13_veg_cs%core%leaf_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%leaf_storage(p)
                  c14_veg_cs%core%leaf_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%core%leaf_storage(p)
               veg_ps%core%leaf_storage(p) = 0._r8
            end if

            ! leaf transfer C and N
            if (abs(veg_cs%core%leaf_xfer(p)) < ccrit) then
               pc = pc + veg_cs%core%leaf_xfer(p)
               veg_cs%core%leaf_xfer(p) = 0._r8
               pn = pn + veg_ns%core%leaf_xfer(p)
               veg_ns%core%leaf_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%leaf_xfer(p)
                  c13_veg_cs%core%leaf_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%leaf_xfer(p)
                  c14_veg_cs%core%leaf_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%core%leaf_xfer(p)
               veg_ps%core%leaf_xfer(p) = 0._r8
            end if

            ! froot C and N
            if (abs(veg_cs%core%froot(p)) < ccrit) then
               pc = pc + veg_cs%core%froot(p)
               veg_cs%core%froot(p) = 0._r8
               pn = pn + veg_ns%core%froot(p)
               veg_ns%core%froot(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%froot(p)
                  c13_veg_cs%core%froot(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%froot(p)
                  c14_veg_cs%core%froot(p) = 0._r8
               endif

               pp = pp + veg_ps%core%froot(p)
               veg_ps%core%froot(p) = 0._r8
            end if

            ! froot storage C and N
            if (abs(veg_cs%core%froot_storage(p)) < ccrit) then
               pc = pc + veg_cs%core%froot_storage(p)
               veg_cs%core%froot_storage(p) = 0._r8
               pn = pn + veg_ns%core%froot_storage(p)
               veg_ns%core%froot_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%froot_storage(p)
                  c13_veg_cs%core%froot_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%froot_storage(p)
                  c14_veg_cs%core%froot_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%core%froot_storage(p)
               veg_ps%core%froot_storage(p) = 0._r8
            end if

            ! froot transfer C and N
            if (abs(veg_cs%core%froot_xfer(p)) < ccrit) then
               pc = pc + veg_cs%core%froot_xfer(p)
               veg_cs%core%froot_xfer(p) = 0._r8
               pn = pn + veg_ns%core%froot_xfer(p)
               veg_ns%core%froot_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%froot_xfer(p)
                  c13_veg_cs%core%froot_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%froot_xfer(p)
                  c14_veg_cs%core%froot_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%core%froot_xfer(p)
               veg_ps%core%froot_xfer(p) = 0._r8
            end if

            if ( crop_prog .and. iscft(veg_pp%itype(p)))then
               ! grain C and N
               if (abs(veg_cs%core%grain(p)) < ccrit) then
                  pc = pc + veg_cs%core%grain(p)
                  veg_cs%core%grain(p) = 0._r8
                  pn = pn + veg_ns%core%grain(p)
                  veg_ns%core%grain(p) = 0._r8
                  pp = pp + veg_ps%core%grain(p)
                  veg_ps%core%grain(p) = 0._r8
               end if

               ! grain storage C and N
               if (abs(veg_cs%core%grain_storage(p)) < ccrit) then
                  pc = pc + veg_cs%core%grain_storage(p)
                  veg_cs%core%grain_storage(p) = 0._r8
                  pn = pn + veg_ns%core%grain_storage(p)
                  veg_ns%core%grain_storage(p) = 0._r8
                  pp = pp + veg_ps%core%grain_storage(p)
                  veg_ps%core%grain_storage(p) = 0._r8
               end if

               ! grain transfer C and N
               if (abs(veg_cs%core%grain_xfer(p)) < ccrit) then
                  pc = pc + veg_cs%core%grain_xfer(p)
                  veg_cs%core%grain_xfer(p) = 0._r8
                  pn = pn + veg_ns%core%grain_xfer(p)
                  veg_ns%core%grain_xfer(p) = 0._r8
                  pp = pp + veg_ps%core%grain_xfer(p)
                  veg_ps%core%grain_xfer(p) = 0._r8
               end if
            end if

            ! livestem C and N
            if (abs(veg_cs%core%livestem(p)) < ccrit) then
               pc = pc + veg_cs%core%livestem(p)
               veg_cs%core%livestem(p) = 0._r8
               pn = pn + veg_ns%core%livestem(p)
               veg_ns%core%livestem(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%livestem(p)
                  c13_veg_cs%core%livestem(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%livestem(p)
                  c14_veg_cs%core%livestem(p) = 0._r8
               endif

               pp = pp + veg_ps%core%livestem(p)
               veg_ps%core%livestem(p) = 0._r8
            end if

            ! livestem storage C and N
            if (abs(veg_cs%core%livestem_storage(p)) < ccrit) then
               pc = pc + veg_cs%core%livestem_storage(p)
               veg_cs%core%livestem_storage(p) = 0._r8
               pn = pn + veg_ns%core%livestem_storage(p)
               veg_ns%core%livestem_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%livestem_storage(p)
                  c13_veg_cs%core%livestem_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%livestem_storage(p)
                  c14_veg_cs%core%livestem_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%core%livestem_storage(p)
               veg_ps%core%livestem_storage(p) = 0._r8
            end if

            ! livestem transfer C and N
            if (abs(veg_cs%core%livestem_xfer(p)) < ccrit) then
               pc = pc + veg_cs%core%livestem_xfer(p)
               veg_cs%core%livestem_xfer(p) = 0._r8
               pn = pn + veg_ns%core%livestem_xfer(p)
               veg_ns%core%livestem_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%livestem_xfer(p)
                  c13_veg_cs%core%livestem_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%livestem_xfer(p)
                  c14_veg_cs%core%livestem_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%core%livestem_xfer(p)
               veg_ps%core%livestem_xfer(p) = 0._r8
            end if

            ! deadstem C and N
            if (abs(veg_cs%core%deadstem(p)) < ccrit) then
               pc = pc + veg_cs%core%deadstem(p)
               veg_cs%core%deadstem(p) = 0._r8
               pn = pn + veg_ns%core%deadstem(p)
               veg_ns%core%deadstem(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%deadstem(p)
                  c13_veg_cs%core%deadstem(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%deadstem(p)
                  c14_veg_cs%core%deadstem(p) = 0._r8
               endif

               pp = pp + veg_ps%core%deadstem(p)
               veg_ps%core%deadstem(p) = 0._r8
            end if

            ! deadstem storage C and N
            if (abs(veg_cs%core%deadstem_storage(p)) < ccrit) then
               pc = pc + veg_cs%core%deadstem_storage(p)
               veg_cs%core%deadstem_storage(p) = 0._r8
               pn = pn + veg_ns%core%deadstem_storage(p)
               veg_ns%core%deadstem_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%deadstem_storage(p)
                  c13_veg_cs%core%deadstem_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%deadstem_storage(p)
                  c14_veg_cs%core%deadstem_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%core%deadstem_storage(p)
               veg_ps%core%deadstem_storage(p) = 0._r8
            end if

            ! deadstem transfer C and N
            if (abs(veg_cs%core%deadstem_xfer(p)) < ccrit) then
               pc = pc + veg_cs%core%deadstem_xfer(p)
               veg_cs%core%deadstem_xfer(p) = 0._r8
               pn = pn + veg_ns%core%deadstem_xfer(p)
               veg_ns%core%deadstem_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%deadstem_xfer(p)
                  c13_veg_cs%core%deadstem_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%deadstem_xfer(p)
                  c14_veg_cs%core%deadstem_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%core%deadstem_xfer(p)
               veg_ps%core%deadstem_xfer(p) = 0._r8
            end if

            ! livecroot C and N
            if (abs(veg_cs%core%livecroot(p)) < ccrit) then
               pc = pc + veg_cs%core%livecroot(p)
               veg_cs%core%livecroot(p) = 0._r8
               pn = pn + veg_ns%core%livecroot(p)
               veg_ns%core%livecroot(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%livecroot(p)
                  c13_veg_cs%core%livecroot(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%livecroot(p)
                  c14_veg_cs%core%livecroot(p) = 0._r8
               endif

               pp = pp + veg_ps%core%livecroot(p)
               veg_ps%core%livecroot(p) = 0._r8
            end if

            ! livecroot storage C and N
            if (abs(veg_cs%core%livecroot_storage(p)) < ccrit) then
               pc = pc + veg_cs%core%livecroot_storage(p)
               veg_cs%core%livecroot_storage(p) = 0._r8
               pn = pn + veg_ns%core%livecroot_storage(p)
               veg_ns%core%livecroot_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%livecroot_storage(p)
                  c13_veg_cs%core%livecroot_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%livecroot_storage(p)
                  c14_veg_cs%core%livecroot_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%core%livecroot_storage(p)
               veg_ps%core%livecroot_storage(p) = 0._r8
            end if

            ! livecroot transfer C and N
            if (abs(veg_cs%core%livecroot_xfer(p)) < ccrit) then
               pc = pc + veg_cs%core%livecroot_xfer(p)
               veg_cs%core%livecroot_xfer(p) = 0._r8
               pn = pn + veg_ns%core%livecroot_xfer(p)
               veg_ns%core%livecroot_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%livecroot_xfer(p)
                  c13_veg_cs%core%livecroot_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%livecroot_xfer(p)
                  c14_veg_cs%core%livecroot_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%core%livecroot_xfer(p)
               veg_ps%core%livecroot_xfer(p) = 0._r8
            end if

            ! deadcroot C and N
            if (abs(veg_cs%core%deadcroot(p)) < ccrit) then
               pc = pc + veg_cs%core%deadcroot(p)
               veg_cs%core%deadcroot(p) = 0._r8
               pn = pn + veg_ns%core%deadcroot(p)
               veg_ns%core%deadcroot(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%deadcroot(p)
                  c13_veg_cs%core%deadcroot(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%deadcroot(p)
                  c14_veg_cs%core%deadcroot(p) = 0._r8
               endif

               pp = pp + veg_ps%core%deadcroot(p)
               veg_ps%core%deadcroot(p) = 0._r8
            end if

            ! deadcroot storage C and N
            if (abs(veg_cs%core%deadcroot_storage(p)) < ccrit) then
               pc = pc + veg_cs%core%deadcroot_storage(p)
               veg_cs%core%deadcroot_storage(p) = 0._r8
               pn = pn + veg_ns%core%deadcroot_storage(p)
               veg_ns%core%deadcroot_storage(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%deadcroot_storage(p)
                  c13_veg_cs%core%deadcroot_storage(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%deadcroot_storage(p)
                  c14_veg_cs%core%deadcroot_storage(p) = 0._r8
               endif

               pp = pp + veg_ps%core%deadcroot_storage(p)
               veg_ps%core%deadcroot_storage(p) = 0._r8
            end if

            ! deadcroot transfer C and N
            if (abs(veg_cs%core%deadcroot_xfer(p)) < ccrit) then
               pc = pc + veg_cs%core%deadcroot_xfer(p)
               veg_cs%core%deadcroot_xfer(p) = 0._r8
               pn = pn + veg_ns%core%deadcroot_xfer(p)
               veg_ns%core%deadcroot_xfer(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%deadcroot_xfer(p)
                  c13_veg_cs%core%deadcroot_xfer(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%deadcroot_xfer(p)
                  c14_veg_cs%core%deadcroot_xfer(p) = 0._r8
               endif

               pp = pp + veg_ps%core%deadcroot_xfer(p)
               veg_ps%core%deadcroot_xfer(p) = 0._r8
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
            if (abs(veg_cs%core%pool(p)) < ccrit) then
               pc = pc + veg_cs%core%pool(p)
               veg_cs%core%pool(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13_veg_cs%core%pool(p)
                  c13_veg_cs%core%pool(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14_veg_cs%core%pool(p)
                  c14_veg_cs%core%pool(p) = 0._r8
               endif
            end if

            if ( crop_prog .and. iscft(veg_pp%itype(p)))then
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
            if (abs(veg_ns%core%pool(p)) < ncrit) then
               pn = pn + veg_ns%core%pool(p)
               veg_ns%core%pool(p) = 0._r8
            end if

            ! ppool (P only)
            if (abs(veg_ps%core%pool(p)) < pcrit) then
               pp = pp + veg_ps%core%pool(p)
               veg_ps%core%pool(p) = 0._r8
            end if

            veg_cs%core%trunc_sink(p) = veg_cs%core%trunc_sink(p) + pc
            veg_ns%core%trunc_sink(p) = veg_ns%core%trunc_sink(p) + pn
            veg_ps%core%trunc_sink(p) = veg_ps%core%trunc_sink(p) + pp

            if ( use_c13 ) then
               c13_veg_cs%core%trunc_sink(p) = c13_veg_cs%core%trunc_sink(p) + pc13
            endif
            if ( use_c14 ) then
               c14_veg_cs%core%trunc_sink(p) = c14_veg_cs%core%trunc_sink(p) + pc14
            endif

         end do ! end of pft loop
      end if ! end of if(.not.use_fates)


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

                  if (abs(col_cs%core%decomp_pools_vr(c,j,k)) < ccrit) then
                     cc = cc + col_cs%core%decomp_pools_vr(c,j,k)
                     col_cs%core%decomp_pools_vr(c,j,k) = 0._r8
                     cn = cn + col_ns%core%decomp_pools_vr(c,j,k)
                     col_ns%core%decomp_pools_vr(c,j,k) = 0._r8
                     if ( use_c13 ) then
                        cc13 = cc13 + c13_col_cs%core%decomp_pools_vr(c,j,k)
                        c13_col_cs%core%decomp_pools_vr(c,j,k) = 0._r8
                     endif
                     if ( use_c14 ) then
                        cc14 = cc14 + c14_col_cs%core%decomp_pools_vr(c,j,k)
                        c14_col_cs%core%decomp_pools_vr(c,j,k) = 0._r8
                     endif
                  end if

               end do

               ! not doing precision control on soil mineral N, since it will
               ! be getting the N truncation flux anyway.

               col_cs%core%trunc_sink_vr(c,j) = col_cs%core%trunc_sink_vr(c,j) + cc
               col_ns%core%trunc_sink_vr(c,j) = col_ns%core%trunc_sink_vr(c,j) + cn
               if ( use_c13 ) then
                  c13_col_cs%core%trunc_sink_vr(c,j) = c13_col_cs%core%trunc_sink_vr(c,j) + cc13
               endif
               if ( use_c14 ) then
                  c14_col_cs%core%trunc_sink_vr(c,j) = c14_col_cs%core%trunc_sink_vr(c,j) + cc14
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
                     if ( col_cs%core%decomp_pools_vr(c,j,l) > 0.0_r8 ) then
                          if(abs(col_cs%core%decomp_pools_vr(c,j,l) / col_ns%core%decomp_pools_vr(c,j,l) - initial_cn_ratio(l) ) > 1.0e-3_r8 &
                          .and. (.not. floating_cn_ratio_decomp_pools(l)) ) then
                        cn_eca = cn_eca - ( col_cs%core%decomp_pools_vr(c,j,l) / initial_cn_ratio(l) - col_ns%core%decomp_pools_vr(c,j,l) )

                          col_ns%core%decomp_pools_vr(c,j,l) = col_cs%core%decomp_pools_vr(c,j,l) / initial_cn_ratio(l)
                     end if
                   end if
                  end do
                  col_ns%core%trunc_sink_vr(c,j) = col_ns%core%trunc_sink_vr(c,j) + cn_eca
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
                     if ( col_ns%core%decomp_pools_vr(c,j,l) < 0.0_r8 .and. floating_cn_ratio_decomp_pools(l) ) then
                        if ( abs(col_ns%core%decomp_pools_vr(c,j,l))  < ncrit ) then
                           cn_eca = cn_eca - ncrit + col_ns%core%decomp_pools_vr(c,j,l)
                           col_ns%core%decomp_pools_vr(c,j,l) = ncrit
                        else
#ifndef _OPENACC                                
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_npools is negative: ',j,l,col_ns%core%decomp_pools_vr(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
                        end if
                     end if
                     if ( col_ps%core%decomp_pools_vr(c,j,l)  < 0.0_r8 .and. floating_cp_ratio_decomp_pools(l) ) then
                        if ( abs(col_ps%core%decomp_pools_vr(c,j,l))  < ncrit/1e4_r8 ) then
                           cp_eca = cp_eca - ncrit/1e4_r8 + col_ps%core%decomp_pools_vr(c,j,l)
                           col_ps%core%decomp_pools_vr(c,j,l) = ncrit/1e4_r8
                         else
#ifndef _OPENACC
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_ppools is negative: ',j,l,col_ps%core%decomp_pools_vr(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
                         end if
                     end if

                  end do

                  col_ns%core%trunc_sink_vr(c,j) = col_ns%core%trunc_sink_vr(c,j) + cn_eca
                  col_ps%core%trunc_sink_vr(c,j) = col_ps%core%trunc_sink_vr(c,j) + cp_eca

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
               if (veg_ns%core%pool(p) < 0._r8) then
                  write(iulog, *) 'error npool_patch is negative: ',p
                  write(iulog, *) 'npool_patch: ', veg_ns%core%pool(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (veg_ps%retransp(p) < 0._r8) then
                  write(iulog, *) 'error retransp_patch is negative: ',p
                  write(iulog, *) 'retransp_patch: ', veg_ps%retransp(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (veg_ps%core%pool(p) < 0._r8) then
                  write(iulog, *) 'error ppool_patch is negative: ',p
                  write(iulog, *) 'ppool_patch: ', veg_ps%core%pool(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
#endif
            end do
         endif

       endif  !if ECA



    end associate

 end subroutine PrecisionControl

end module PrecisionControlMod
