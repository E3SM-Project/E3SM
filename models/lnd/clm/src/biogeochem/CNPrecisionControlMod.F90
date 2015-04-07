module CNPrecisionControlMod

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables 
  ! 
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use CNVegCarbonStateType   , only : cnveg_carbonstate_type
  use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
  use PatchType              , only : patch
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNPrecisionControl
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine CNPrecisionControl(num_soilp, filter_soilp, &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       cnveg_nitrogenstate_inst)
    !
    ! !DESCRIPTION: 
    ! Force leaf and deadstem c and n to 0 if they get too small.
    !
    ! !USES:
    use clm_varctl , only : iulog, use_c13, use_c14
    use clm_varpar , only : crop_prog
    use pftconMod  , only : nc3crop
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p,j,k  ! indices
    integer :: fp       ! filter indices
    real(r8):: pc,pn    ! truncation terms for patch-level corrections
    real(r8):: pc13     ! truncation terms for patch-level corrections
    real(r8):: pc14     ! truncation terms for patch-level corrections
    real(r8):: cc14     ! truncation terms for column-level corrections
    real(r8):: ccrit    ! critical carbon state value for truncation
    real(r8):: ncrit    ! critical nitrogen state value for truncation
    !-----------------------------------------------------------------------

    ! cnveg_carbonstate_inst%cpool_patch                     Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! cnveg_carbonstate_inst%deadcrootc_patch                Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! cnveg_carbonstate_inst%deadcrootc_storage_patch        Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! cnveg_carbonstate_inst%deadcrootc_xfer_patch           Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! cnveg_carbonstate_inst%deadstemc_patch                 Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! cnveg_carbonstate_inst%deadstemc_storage_patch         Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! cnveg_carbonstate_inst%deadstemc_xfer_patch            Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! cnveg_carbonstate_inst%frootc_patch                    Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! cnveg_carbonstate_inst%frootc_storage_patch            Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! cnveg_carbonstate_inst%frootc_xfer_patch               Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! cnveg_carbonstate_inst%gresp_storage_patch             Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! cnveg_carbonstate_inst%gresp_xfer_patch                Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! cnveg_carbonstate_inst%leafc_patch                     Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! cnveg_carbonstate_inst%leafc_storage_patch             Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! cnveg_carbonstate_inst%leafc_xfer_patch                Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! cnveg_carbonstate_inst%livecrootc_patch                Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! cnveg_carbonstate_inst%livecrootc_storage_patch        Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! cnveg_carbonstate_inst%livecrootc_xfer_patch           Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! cnveg_carbonstate_inst%livestemc_patch                 Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! cnveg_carbonstate_inst%livestemc_storage_patch         Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! cnveg_carbonstate_inst%livestemc_xfer_patch            Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! cnveg_carbonstate_inst%ctrunc_patch                    Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    ! cnveg_carbonstate_inst%xsmrpool_patch                  Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool                  
    ! cnveg_carbonstate_inst%grainc_patch                    Output:  [real(r8) (:)     ]  (gC/m2) grain C                                   
    ! cnveg_carbonstate_inst%grainc_storage_patch            Output:  [real(r8) (:)     ]  (gC/m2) grain C storage                           
    ! cnveg_carbonstate_inst%grainc_xfer_patch               Output:  [real(r8) (:)     ]  (gC/m2) grain C transfer                          
    
    ! cnveg_nitrogenstate_inst%deadcrootn_patch              Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
    ! cnveg_nitrogenstate_inst%deadcrootn_storage_patch      Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
    ! cnveg_nitrogenstate_inst%deadcrootn_xfer_patch         Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
    ! cnveg_nitrogenstate_inst%deadstemn_patch               Output:  [real(r8) (:)     ]  (gN/m2) dead stem N                               
    ! cnveg_nitrogenstate_inst%deadstemn_storage_patch       Output:  [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
    ! cnveg_nitrogenstate_inst%deadstemn_xfer_patch          Output:  [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
    ! cnveg_nitrogenstate_inst%frootn_patch                  Output:  [real(r8) (:)     ]  (gN/m2) fine root N                               
    ! cnveg_nitrogenstate_inst%frootn_storage_patch          Output:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
    ! cnveg_nitrogenstate_inst%frootn_xfer_patch             Output:  [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
    ! cnveg_nitrogenstate_inst%leafn_patch                   Output:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
    ! cnveg_nitrogenstate_inst%leafn_storage_patch           Output:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
    ! cnveg_nitrogenstate_inst%leafn_xfer_patch              Output:  [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
    ! cnveg_nitrogenstate_inst%livecrootn_patch              Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
    ! cnveg_nitrogenstate_inst%livecrootn_storage_patch      Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
    ! cnveg_nitrogenstate_inst%livecrootn_xfer_patch         Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
    ! cnveg_nitrogenstate_inst%grainn_patch                  Output:  [real(r8) (:)     ]  (gC/m2) grain N                                   
    ! cnveg_nitrogenstate_inst%grainn_storage_patch          Output:  [real(r8) (:)     ]  (gC/m2) grain N storage                           
    ! cnveg_nitrogenstate_inst%grainn_xfer_patch             Output:  [real(r8) (:)     ]  (gC/m2) grain N transfer                          
    ! cnveg_nitrogenstate_inst%livestemn_patch               Output:  [real(r8) (:)     ]  (gN/m2) live stem N                               
    ! cnveg_nitrogenstate_inst%livestemn_storage_patch       Output:  [real(r8) (:)     ]  (gN/m2) live stem N storage                       
    ! cnveg_nitrogenstate_inst%livestemn_xfer_patch          Output:  [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
    ! cnveg_nitrogenstate_inst%npool_patch                   Output:  [real(r8) (:)     ]  (gN/m2) temporary plant N pool                    
    ! cnveg_nitrogenstate_inst%ntrunc_patch                  Output:  [real(r8) (:)     ]  (gN/m2) patch-level sink for N truncation           
    ! cnveg_nitrogenstate_inst%retransn_patch                Output:  [real(r8) (:)     ]  (gN/m2) plant pool of retranslocated N            
    
    associate(                                           &
         cs     => cnveg_carbonstate_inst              , &
         ns     => cnveg_nitrogenstate_inst            , &
         c13cs  => c13_cnveg_carbonstate_inst          , &
         c14cs  => c14_cnveg_carbonstate_inst            &
         )

      ! set the critical carbon state value for truncation (gC/m2)
      ccrit = 1.e-8_r8

      ! set the critical nitrogen state value for truncation (gN/m2)
      ncrit = 1.e-8_r8

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! initialize the patch-level C and N truncation terms
         pc = 0._r8
         pn = 0._r8
         if ( use_c13 ) pc13 = 0._r8
         if ( use_c14 ) pc14 = 0._r8

         ! do tests on state variables for precision control
         ! for linked C-N state variables, perform precision test on
         ! the C component, but truncate C, C13, and N components

         ! leaf C and N
         if (abs(cs%leafc_patch(p)) < ccrit) then
            pc = pc + cs%leafc_patch(p)
            cs%leafc_patch(p) = 0._r8

            pn = pn + ns%leafn_patch(p)
            ns%leafn_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%leafc_patch(p)
               c13cs%leafc_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%leafc_patch(p)
               c14cs%leafc_patch(p) = 0._r8
            endif
         end if

         ! leaf storage C and N
         if (abs(cs%leafc_storage_patch(p)) < ccrit) then
            pc = pc + cs%leafc_storage_patch(p)
            cs%leafc_storage_patch(p) = 0._r8

            pn = pn + ns%leafn_storage_patch(p)
            ns%leafn_storage_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%leafc_storage_patch(p)
               c13cs%leafc_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%leafc_storage_patch(p)
               c14cs%leafc_storage_patch(p) = 0._r8
            endif
         end if

         ! leaf transfer C and N
         if (abs(cs%leafc_xfer_patch(p)) < ccrit) then
            pc = pc + cs%leafc_xfer_patch(p)
            cs%leafc_xfer_patch(p) = 0._r8

            pn = pn + ns%leafn_xfer_patch(p)
            ns%leafn_xfer_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%leafc_xfer_patch(p)
               c13cs%leafc_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%leafc_xfer_patch(p)
               c14cs%leafc_xfer_patch(p) = 0._r8
            endif
         end if

         ! froot C and N
         if (abs(cs%frootc_patch(p)) < ccrit) then
            pc = pc + cs%frootc_patch(p)
            cs%frootc_patch(p) = 0._r8

            pn = pn + ns%frootn_patch(p)
            ns%frootn_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%frootc_patch(p)
               c13cs%frootc_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%frootc_patch(p)
               c14cs%frootc_patch(p) = 0._r8
            endif
         end if

         ! froot storage C and N
         if (abs(cs%frootc_storage_patch(p)) < ccrit) then
            pc = pc + cs%frootc_storage_patch(p)
            cs%frootc_storage_patch(p) = 0._r8

            pn = pn + ns%frootn_storage_patch(p)
            ns%frootn_storage_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%frootc_storage_patch(p)
               c13cs%frootc_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%frootc_storage_patch(p)
               c14cs%frootc_storage_patch(p) = 0._r8
            endif
         end if

         ! froot transfer C and N
         if (abs(cs%frootc_xfer_patch(p)) < ccrit) then
            pc = pc + cs%frootc_xfer_patch(p)
            cs%frootc_xfer_patch(p) = 0._r8

            pn = pn + ns%frootn_xfer_patch(p)
            ns%frootn_xfer_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%frootc_xfer_patch(p)
               c13cs%frootc_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%frootc_xfer_patch(p)
               c14cs%frootc_xfer_patch(p) = 0._r8
            endif
         end if

         if ( crop_prog .and. patch%itype(p) >= nc3crop )then
            ! grain C and N
            if (abs(cs%grainc_patch(p)) < ccrit) then
               pc = pc + cs%grainc_patch(p)
               cs%grainc_patch(p) = 0._r8

               pn = pn + ns%grainn_patch(p)
               ns%grainn_patch(p) = 0._r8
            end if

            ! grain storage C and N
            if (abs(cs%grainc_storage_patch(p)) < ccrit) then
               pc = pc + cs%grainc_storage_patch(p)
               cs%grainc_storage_patch(p) = 0._r8

               pn = pn + ns%grainn_storage_patch(p)
               ns%grainn_storage_patch(p) = 0._r8
            end if

            ! grain transfer C and N
            if (abs(cs%grainc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%grainc_xfer_patch(p)
               cs%grainc_xfer_patch(p) = 0._r8

               pn = pn + ns%grainn_xfer_patch(p)
               ns%grainn_xfer_patch(p) = 0._r8
            end if
         end if

         ! livestem C and N
         if (abs(cs%livestemc_patch(p)) < ccrit) then
            pc = pc + cs%livestemc_patch(p)
            cs%livestemc_patch(p) = 0._r8

            pn = pn + ns%livestemn_patch(p)
            ns%livestemn_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%livestemc_patch(p)
               c13cs%livestemc_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%livestemc_patch(p)
               c14cs%livestemc_patch(p) = 0._r8
            endif
         end if

         ! livestem storage C and N
         if (abs(cs%livestemc_storage_patch(p)) < ccrit) then
            pc = pc + cs%livestemc_storage_patch(p)
            cs%livestemc_storage_patch(p) = 0._r8

            pn = pn + ns%livestemn_storage_patch(p)
            ns%livestemn_storage_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%livestemc_storage_patch(p)
               c13cs%livestemc_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%livestemc_storage_patch(p)
               c14cs%livestemc_storage_patch(p) = 0._r8
            endif
         end if

         ! livestem transfer C and N
         if (abs(cs%livestemc_xfer_patch(p)) < ccrit) then
            pc = pc + cs%livestemc_xfer_patch(p)
            cs%livestemc_xfer_patch(p) = 0._r8

            pn = pn + ns%livestemn_xfer_patch(p)
            ns%livestemn_xfer_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%livestemc_xfer_patch(p)
               c13cs%livestemc_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%livestemc_xfer_patch(p)
               c14cs%livestemc_xfer_patch(p) = 0._r8
            endif
         end if

         ! deadstem C and N
         if (abs(cs%deadstemc_patch(p)) < ccrit) then
            pc = pc + cs%deadstemc_patch(p)
            cs%deadstemc_patch(p) = 0._r8

            pn = pn + ns%deadstemn_patch(p)
            ns%deadstemn_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%deadstemc_patch(p)
               c13cs%deadstemc_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%deadstemc_patch(p)
               c14cs%deadstemc_patch(p) = 0._r8
            endif
         end if

         ! deadstem storage C and N
         if (abs(cs%deadstemc_storage_patch(p)) < ccrit) then
            pc = pc + cs%deadstemc_storage_patch(p)
            cs%deadstemc_storage_patch(p) = 0._r8

            pn = pn + ns%deadstemn_storage_patch(p)
            ns%deadstemn_storage_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%deadstemc_storage_patch(p)
               c13cs%deadstemc_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%deadstemc_storage_patch(p)
               c14cs%deadstemc_storage_patch(p) = 0._r8
            endif
         end if

         ! deadstem transfer C and N
         if (abs(cs%deadstemc_xfer_patch(p)) < ccrit) then
            pc = pc + cs%deadstemc_xfer_patch(p)
            cs%deadstemc_xfer_patch(p) = 0._r8

            pn = pn + ns%deadstemn_xfer_patch(p)
            ns%deadstemn_xfer_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%deadstemc_xfer_patch(p)
               c13cs%deadstemc_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%deadstemc_xfer_patch(p)
               c14cs%deadstemc_xfer_patch(p) = 0._r8
            endif
         end if

         ! livecroot C and N
         if (abs(cs%livecrootc_patch(p)) < ccrit) then
            pc = pc + cs%livecrootc_patch(p)
            cs%livecrootc_patch(p) = 0._r8

            pn = pn + ns%livecrootn_patch(p)
            ns%livecrootn_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%livecrootc_patch(p)
               c13cs%livecrootc_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%livecrootc_patch(p)
               c14cs%livecrootc_patch(p) = 0._r8
            endif
         end if

         ! livecroot storage C and N
         if (abs(cs%livecrootc_storage_patch(p)) < ccrit) then
            pc = pc + cs%livecrootc_storage_patch(p)
            cs%livecrootc_storage_patch(p) = 0._r8

            pn = pn + ns%livecrootn_storage_patch(p)
            ns%livecrootn_storage_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%livecrootc_storage_patch(p)
               c13cs%livecrootc_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%livecrootc_storage_patch(p)
               c14cs%livecrootc_storage_patch(p) = 0._r8
            endif
         end if

         ! livecroot transfer C and N
         if (abs(cs%livecrootc_xfer_patch(p)) < ccrit) then
            pc = pc + cs%livecrootc_xfer_patch(p)
            cs%livecrootc_xfer_patch(p) = 0._r8

            pn = pn + ns%livecrootn_xfer_patch(p)
            ns%livecrootn_xfer_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%livecrootc_xfer_patch(p)
               c13cs%livecrootc_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%livecrootc_xfer_patch(p)
               c14cs%livecrootc_xfer_patch(p) = 0._r8
            endif
         end if

         ! deadcroot C and N
         if (abs(cs%deadcrootc_patch(p)) < ccrit) then
            pc = pc + cs%deadcrootc_patch(p)
            cs%deadcrootc_patch(p) = 0._r8

            pn = pn + ns%deadcrootn_patch(p)
            ns%deadcrootn_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%deadcrootc_patch(p)
               c13cs%deadcrootc_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%deadcrootc_patch(p)
               c14cs%deadcrootc_patch(p) = 0._r8
            endif
         end if

         ! deadcroot storage C and N
         if (abs(cs%deadcrootc_storage_patch(p)) < ccrit) then
            pc = pc + cs%deadcrootc_storage_patch(p)
            cs%deadcrootc_storage_patch(p) = 0._r8

            pn = pn + ns%deadcrootn_storage_patch(p)
            ns%deadcrootn_storage_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%deadcrootc_storage_patch(p)
               c13cs%deadcrootc_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%deadcrootc_storage_patch(p)
               c14cs%deadcrootc_storage_patch(p) = 0._r8
            endif
         end if

         ! deadcroot transfer C and N
         if (abs(cs%deadcrootc_xfer_patch(p)) < ccrit) then
            pc = pc + cs%deadcrootc_xfer_patch(p)
            cs%deadcrootc_xfer_patch(p) = 0._r8

            pn = pn + ns%deadcrootn_xfer_patch(p)
            ns%deadcrootn_xfer_patch(p) = 0._r8

            if ( use_c13 ) then
               pc13 = pc13 + c13cs%deadcrootc_xfer_patch(p)
               c13cs%deadcrootc_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%deadcrootc_xfer_patch(p)
               c14cs%deadcrootc_xfer_patch(p) = 0._r8
            endif
         end if

         ! gresp_storage (C only)
         if (abs(cs%gresp_storage_patch(p)) < ccrit) then
            pc = pc + cs%gresp_storage_patch(p)
            cs%gresp_storage_patch(p) = 0._r8
            if ( use_c13 ) then
               pc13 = pc13 + c13cs%gresp_storage_patch(p)
               c13cs%gresp_storage_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%gresp_storage_patch(p)
               c14cs%gresp_storage_patch(p) = 0._r8
            endif
         end if

         ! gresp_xfer(c only)
         if (abs(cs%gresp_xfer_patch(p)) < ccrit) then
            pc = pc + cs%gresp_xfer_patch(p)
            cs%gresp_xfer_patch(p) = 0._r8
            if ( use_c13 ) then
               pc13 = pc13 + c13cs%gresp_xfer_patch(p)
               c13cs%gresp_xfer_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%gresp_xfer_patch(p)
               c14cs%gresp_xfer_patch(p) = 0._r8
            endif
         end if

         ! cpool (C only)
         if (abs(cs%cpool_patch(p)) < ccrit) then
            pc = pc + cs%cpool_patch(p)
            cs%cpool_patch(p) = 0._r8
            if ( use_c13 ) then
               pc13 = pc13 + c13cs%cpool_patch(p)
               c13cs%cpool_patch(p) = 0._r8
            endif
            if ( use_c14 ) then
               pc14 = pc14 + c14cs%cpool_patch(p)
               c14cs%cpool_patch(p) = 0._r8
            endif
         end if

         if ( crop_prog .and. patch%itype(p) >= nc3crop )then
            ! xsmrpool (C only)
            if (abs(cs%xsmrpool_patch(p)) < ccrit) then
               pc = pc + cs%xsmrpool_patch(p)
               cs%xsmrpool_patch(p) = 0._r8
            end if
         end if

         ! retransn (N only)
         if (abs(ns%retransn_patch(p)) < ncrit) then
            pn = pn + ns%retransn_patch(p)
            ns%retransn_patch(p) = 0._r8
         end if

         ! npool (N only)
         if (abs(ns%npool_patch(p)) < ncrit) then
            pn = pn + ns%npool_patch(p)
            ns%npool_patch(p) = 0._r8
         end if

         cs%ctrunc_patch(p) = cs%ctrunc_patch(p) + pc

         ns%ntrunc_patch(p) = ns%ntrunc_patch(p) + pn

         if ( use_c13 ) then
            c13cs%ctrunc_patch(p) = c13cs%ctrunc_patch(p) + pc13
         endif
         if ( use_c14 ) then
            c14cs%ctrunc_patch(p) = c14cs%ctrunc_patch(p) + pc14
         endif

      end do ! end of patch loop

    end associate

 end subroutine CNPrecisionControl

end module CNPrecisionControlMod
