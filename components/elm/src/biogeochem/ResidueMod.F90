module ResidueMod
  !-----------------------------------------------------------------------
  ! !MODULE: ResidueMod.F90
  !
  ! !DESCRIPTION:
  ! Module holding routines for the change of surface residue C/N/P pools.

  !
  ! !USES:
  use elm_varpar              , only : i_met_lit, i_cel_lit, i_lig_lit
  use elm_varpar              , only : nlit_pools, nlevdecomp
  use ColumnDataType          , only : col_cs, col_ns, col_ps
  use ColumnDataType          , only : col_cf, col_nf, col_pf
  use CNStateType             , only : cnstate_type
  use CropType                , only : crop_type
  use VegetationType          , only : veg_pp

  use timeinfoMod
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNResidueToColumn
  !

contains

  !-----------------------------------------------------------------------
  subroutine CNResidueToColumn (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       crop_vars, cnstate_vars)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to assign residue to the three litter pools
    !
    ! !USES:
      !$acc routine seq
    use elm_varcon  , only : krsd
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(crop_type)         , intent(inout) :: crop_vars
    type(cnstate_type)      , intent(in)    :: cnstate_vars

    ! !LOCAL VARIABLES:
    integer  :: fc,fp,c,p       ! indices
    integer  :: k,j
    integer  :: jday
    real(r8) :: litr_prof(1:nlevdecomp)
    real(r8) :: wt_col
    real(r8) :: dtime
    !-----------------------------------------------------------------------

    associate(                                                                                       &
       residue_prof                        =>    cnstate_vars%residue_prof_patch,     & ! Input: [real(r8) (:,:) ] (1/m) profile of residue litter input
       tillage_prof                        =>    cnstate_vars%tillage_prof_patch,     & ! Input: [real(r8) (:,:) ] (1/m) profile of residue litter input by tillage

       residue_cpools                      =>    col_cs%residue_cpools,               & ! Input: [real(r8) (:,:) ] (gC/m2) surface residue (surface litter) c pools
       residue_npools                      =>    col_ns%residue_npools,               & ! Input: [real(r8) (:,:) ] (gN/m2) surface residue (surface litter) n pools
       residue_ppools                      =>    col_ps%residue_ppools,               & ! Input: [real(r8) (:,:) ] (gP/m2) surface residue (surface litter) p pools

       tilldate                            =>    crop_vars%tilldate_patch,            & ! Inout: [integer  (:) ]  tillage date
       istilled                            =>    crop_vars%istilled_patch,            & ! Inout: [logical  (:) ]  Flag, true if tillage performed
       tilleffic                           =>    crop_vars%tilleffic_patch,           & ! Input: [real(r8) (:) ]  tillage efficiency

       residue2litr                        =>    cnstate_vars%residue2litr_patch,     & ! Output: [real(r8) (:) ] Residue to litter conversion rate (1/s)
       
       residue_to_litr_met_c               =>    col_cf%residue_to_litr_met_c,        & ! Output: [real(r8) (:,:) ] C fluxes associated with residue to litter metabolic pool (gC/m3/s)
       residue_to_litr_cel_c               =>    col_cf%residue_to_litr_cel_c,        & ! Output: [real(r8) (:,:) ] C fluxes associated with residue to litter cellulose pool (gC/m3/s) 
       residue_to_litr_lig_c               =>    col_cf%residue_to_litr_lig_c,        & ! Output: [real(r8) (:,:) ] C fluxes associated with residue to litter lignin pool (gC/m3/s)
       
       residue_to_litr_met_n               =>    col_nf%residue_to_litr_met_n,        & ! Output: [real(r8) (:,:) ] N fluxes associated with residue to litter metabolic pool (gN/m3/s) 
       residue_to_litr_cel_n               =>    col_nf%residue_to_litr_cel_n,        & ! Output: [real(r8) (:,:) ] N fluxes associated with residue to litter cellulose pool (gN/m3/s)
       residue_to_litr_lig_n               =>    col_nf%residue_to_litr_lig_n,        & ! Output: [real(r8) (:,:) ] N fluxes associated with residue to litter lignin pool (gN/m3/s)

       residue_to_litr_met_p               =>    col_pf%residue_to_litr_met_p,        & ! Output: [real(r8) (:,:) ] P fluxes associated with residue to litter metabolic pool (gP/m3/s)
       residue_to_litr_cel_p               =>    col_pf%residue_to_litr_cel_p,        & ! Output: [real(r8) (:,:) ] P fluxes associated with residue to litter cellulose pool (gP/m3/s)
       residue_to_litr_lig_p               =>    col_pf%residue_to_litr_lig_p         & ! Output: [real(r8) (:,:) ] P fluxes associated with residue to litter lignin pool (gP/m3/s)
       )

    ! Get time step
    dtime = dtime_mod
    jday  = jday_mod

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       
       residue_to_litr_met_c(c,:) = 0._r8
       residue_to_litr_cel_c(c,:) = 0._r8
       residue_to_litr_lig_c(c,:) = 0._r8

       residue_to_litr_met_n(c,:) = 0._r8
       residue_to_litr_cel_n(c,:) = 0._r8
       residue_to_litr_lig_n(c,:) = 0._r8

       residue_to_litr_met_p(c,:) = 0._r8
       residue_to_litr_cel_p(c,:) = 0._r8
       residue_to_litr_lig_p(c,:) = 0._r8
    end do

    do fp = 1, num_soilp
       p = filter_soilp(fp)
       c = veg_pp%column(p)
       wt_col = veg_pp%wtcol(p)
       
       ! update residue to litter conversion rate
       if ( (tilldate(p)==jday) .and. (.not. istilled(p)) ) then
          ! tillage time
          residue2litr(p) = tilleffic(p) * 1._r8/dtime
          litr_prof = tillage_prof(p,1:nlevdecomp) 
          istilled(p) = .True.
       else
          ! no tillage time
          residue2litr(p) = min( krsd, 1._r8/dtime )
          litr_prof = residue_prof(p,1:nlevdecomp)
       end if  
       
       ! residue C/N/P to litter C/N/P
       residue_to_litr_met_c(c,1:nlevdecomp) = residue_to_litr_met_c(c,1:nlevdecomp) + &
          residue_cpools(p,i_met_lit) * litr_prof * residue2litr(p) * wt_col
       residue_to_litr_cel_c(c,1:nlevdecomp) = residue_to_litr_cel_c(c,1:nlevdecomp) + &
          residue_cpools(p,i_cel_lit) * litr_prof * residue2litr(p) * wt_col
       residue_to_litr_lig_c(c,1:nlevdecomp) = residue_to_litr_lig_c(c,1:nlevdecomp) + &
          residue_cpools(p,i_lig_lit) * litr_prof * residue2litr(p) * wt_col

       residue_to_litr_met_n(c,1:nlevdecomp) = residue_to_litr_met_n(c,1:nlevdecomp) + &
          residue_npools(p,i_met_lit) * litr_prof * residue2litr(p) * wt_col
       residue_to_litr_cel_n(c,1:nlevdecomp) = residue_to_litr_cel_n(c,1:nlevdecomp) + &
          residue_npools(p,i_cel_lit) * litr_prof * residue2litr(p) * wt_col
       residue_to_litr_lig_n(c,1:nlevdecomp) = residue_to_litr_lig_n(c,1:nlevdecomp) + &
          residue_npools(p,i_lig_lit) * litr_prof * residue2litr(p) * wt_col

       residue_to_litr_met_p(c,1:nlevdecomp) = residue_to_litr_met_p(c,1:nlevdecomp) + &
          residue_ppools(p,i_met_lit) * litr_prof * residue2litr(p) * wt_col
       residue_to_litr_cel_p(c,1:nlevdecomp) = residue_to_litr_cel_p(c,1:nlevdecomp) + &
          residue_ppools(p,i_cel_lit) * litr_prof * residue2litr(p) * wt_col
       residue_to_litr_lig_p(c,1:nlevdecomp) = residue_to_litr_lig_p(c,1:nlevdecomp) + &
          residue_ppools(p,i_lig_lit) * litr_prof * residue2litr(p) * wt_col
    end do

    end associate

  end subroutine CNResidueToColumn

end module ResidueMod
