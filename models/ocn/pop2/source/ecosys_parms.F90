MODULE ecosys_parms

  !-----------------------------------------------------------------------------
  !   This module manages the parameter variables for the module ecosys_mod.
  !   Most of the variables are not parameters in the Fortran sense. In the
  !   the Fortran sense, they are vanilla module variables.
  !
  !   This modules handles initializing the variables to default values and
  !   reading them from the namelist ecosys_parms. The values used are echoed
  !   to stdout for record keeping purposes.
  !
  !   CVS:$Id: ecosys_parms.F90 941 2006-05-12 21:36:48Z klindsay $
  !   CVS:$Name$
  !-----------------------------------------------------------------------------
  !   Modified to include parameters for diazotrophs, JKM  4/2002
  !-----------------------------------------------------------------------------
  !   variables/subroutines/function used from other modules
  !   The following are used extensively in this ecosys, so are used at
  !   the module level. The use statements for variables that are only needed
  !   locally are located at the module subprogram level.
  !-----------------------------------------------------------------------------

  USE exit_mod, ONLY : sigAbort, exit_POP
  USE communicate, ONLY : my_task, master_task
  USE constants, ONLY : c1
  USE kinds_mod
  USE io_tools, ONLY : document

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !   all module variables are public and should have their values preserved
  !-----------------------------------------------------------------------------

  PUBLIC
  SAVE

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       spd = 86400.0_r8,        & ! number of seconds in a day
       dps = c1 / spd,          & ! number of days in a second
       yps = c1 / (365.0_r8*spd)  ! number of years in a second

  !-----------------------------------------------------------------------------
  !   functional group
  !-----------------------------------------------------------------------------

  TYPE, PUBLIC :: autotroph_type
     CHARACTER(char_len) :: sname, lname
     LOGICAL(KIND=log_kind) :: &
        Nfixer,                             & ! flag set to true if this autotroph fixes N2
        imp_calcifier,                      & ! flag set to true if this autotroph implicitly handles calcification
        exp_calcifier                         ! flag set to true if this autotroph explicitly handles calcification
     INTEGER (KIND=int_kind) :: &
        Chl_ind, C_ind, Fe_ind,             & ! tracer indices for Chl, C, Fe content
        Si_ind, CaCO3_ind,                  & ! tracer indices for Si, CaCO3 content
        grazee_ind                            ! which grazee category does autotroph belong to
     REAL(KIND=r8) :: &
        kFe, kPO4, kDOP, kNO3, kNH4, kSiO3, & ! nutrient uptake half-sat constants
        Qp,                                 & ! P/C ratio
        gQfe_0, gQfe_min,                   & ! initial and minimum fe/C ratio
        alphaPI,                            & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
        PCref,                              & ! max C-spec. grth rate at tref (1/sec)
        thetaN_max,                         & ! max thetaN (Chl/N) (mg Chl/mmol N)
        loss_thres, loss_thres2,            & ! conc. where losses go to zero
        temp_thres,                         & ! Temp. where concentration threshold and photosynth. rate drops
        mort, mort2,                        & ! linear and quadratic mortality rates (1/sec), (1/sec/((mmol C/m3))
        agg_rate_max, agg_rate_min,         & ! max and min agg. rate (1/d)
        z_umax_0,                           & ! max zoo growth rate at tref (1/sec)
        z_grz,                              & ! grazing coef. (mmol C/m^3)
        graze_zoo, graze_poc, graze_doc,    & ! routing of grazed term, remainder goes to dic
        loss_poc,                           & ! routing of loss term
        f_zoo_detr                            ! fraction of zoo losses to detrital 
  END TYPE

  INTEGER (KIND=int_kind), PARAMETER :: &
     autotroph_cnt   = 3, &
     sp_ind          = 1, &  ! small phytoplankton
     diat_ind        = 2, &  ! diatoms
     diaz_ind        = 3     ! diazotrophs

  TYPE(autotroph_type), DIMENSION(autotroph_cnt) :: autotrophs

  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       parm_Red_D_C_P  = 117.0_r8,                 & ! carbon:phosphorus
       parm_Red_D_N_P  =  16.0_r8,                 & ! nitrogen:phosphorus
       parm_Red_D_O2_P = 170.0_r8,                 & ! oxygen:phosphorus
       parm_Remin_D_O2_P = 138.0_r8,               & ! oxygen:phosphorus
       parm_Red_P_C_P  = parm_Red_D_C_P,                 & ! carbon:phosphorus
       parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,  & ! carbon:nitrogen
       parm_Red_P_C_N  = parm_Red_D_C_N,                 & ! carbon:nitrogen
       parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P, & ! carbon:oxygen
       parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P, & ! carbon:oxygen
       parm_Red_P_C_O2 = parm_Red_D_C_O2,                & ! carbon:oxygen
       parm_Red_Fe_C   = 3.0e-6_r8,                & ! iron:carbon
       parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_r8! carbon:oxygen
                                                           ! for diazotrophs

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via namelist input
  !----------------------------------------------------------------------------

  REAL(KIND=r8) :: &
       parm_Fe_bioavail,      & ! fraction of Fe flux that is bioavailable
       parm_o2_min,           & ! min O2 needed for prod & consump. (nmol/cm^3)
       parm_o2_min_delta,     & ! width of min O2 range (nmol/cm^3)
       parm_kappa_nitrif,     & ! nitrification inverse time constant (1/sec)
       parm_nitrif_par_lim,   & ! PAR limit for nitrif. (W/m^2)
       parm_z_mort_0,         & ! zoo linear mort rate (1/sec)
       parm_z_mort2_0,        & ! zoo quad mort rate (1/sec/((mmol C/m3))
       parm_labile_ratio,     & ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
       parm_POMbury,          & ! scale factor for burial of POC, PON, and POP
       parm_BSIbury,          & ! scale factor burial of bSi
       parm_fe_scavenge_rate0,& ! base scavenging rate
       parm_f_prod_sp_CaCO3,  & !fraction of sp prod. as CaCO3 prod.
       parm_POC_diss,         & ! base POC diss len scale
       parm_SiO2_diss,        & ! base SiO2 diss len scale
       parm_CaCO3_diss          ! base CaCO3 diss len scale

  REAL(KIND=r8), DIMENSION(4) :: &
       parm_scalelen_z,       & ! depths of prescribed scalelen values
       parm_scalelen_vals       ! prescribed scalelen values

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       fe_scavenge_thres1 = 0.8e-3_r8,  & !upper thres. for Fe scavenging
       dust_fescav_scale  = 1.0e9,      & !dust scavenging scale factor
       fe_max_scale2      = 1200.0_r8     !unitless scaling coeff.

  !---------------------------------------------------------------------
  !     Compute iron remineralization and flux out.
  !     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
  !                        --------- *  ---------- * ----------
  !			    gDust       55.847 gFe     molFe
  !
  !     dust_to_Fe          conversion - dust to iron (nmol Fe/g Dust) 
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       dust_to_Fe=0.035_r8/55.847_r8*1.0e9_r8
 
  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
      caco3_poc_min    = 0.4_r8,  & !minimum proportionality between 
                                          !   QCaCO3 and grazing losses to POC 
                                          !   (mmol C/mmol CaCO3)
      spc_poc_fac      = 0.11_r8, & !small phyto grazing factor (1/mmolC)
      f_graze_sp_poc_lim = 0.3_r8, & 
      f_photosp_CaCO3  = 0.4_r8,  & !proportionality between small phyto 
                                          !    production and CaCO3 production
      f_graze_CaCO3_remin = 0.33_r8, & !fraction of spCaCO3 grazing 
                                             !          which is remin
      f_graze_si_remin    = 0.35_r8      !fraction of diatom Si grazing 
                                             !          which is remin

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       r_Nfix_photo=1.25_r8         ! N fix relative to C fix (non-dim)

  !-----------------------------------------------------------------------
  !     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
  !     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
  !     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
  !-----------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
      Q             = 0.137_r8,  & !N/C ratio (mmol/mmol) of phyto & zoo
      Qp_zoo_pom    = 0.00855_r8,& !P/C ratio (mmol/mmol) zoo & pom
      Qfe_zoo       = 3.0e-6_r8, & !zooplankton fe/C ratio
      gQsi_0        = 0.137_r8,  & !initial Si/C ratio
      gQsi_max      = 0.685_r8,  & !max Si/C ratio
      gQsi_min      = 0.0457_r8, & !min Si/C ratio
      QCaCO3_max    = 0.4_r8,    & !max QCaCO3
      ! carbon:nitrogen ratio for denitrification
      denitrif_C_N  = parm_Red_D_C_P/136.0_r8

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
      thres_z1          = 100.0e2_r8, & !threshold = C_loss_thres for z shallower than this (cm)
      thres_z2          = 150.0e2_r8, & !threshold = 0 for z deeper than this (cm)
      loss_thres_zoo    = 0.005_r8,  & !zoo conc. where losses go to zero
      CaCO3_temp_thres1 = 6.0_r8,   & !upper temp threshold for CaCO3 prod
      CaCO3_temp_thres2 = -2.0_r8,  & !lower temp threshold
      CaCO3_sp_thres    = 4.0_r8      ! bloom condition thres (mmolC/m3)

  !---------------------------------------------------------------------
  !     fraction of incoming shortwave assumed to be PAR
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       f_qsw_par = 0.45_r8   ! PAR fraction

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       Tref = 30.0_r8, & ! reference temperature (C)
       Q_10 = 1.5_r8     ! factor for temperature dependence (non-dim)

  !---------------------------------------------------------------------
  !  DOM parameters for refractory components and DOP uptake
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       DOC_reminR  = (c1/250.0_r8) * dps,          & ! rate for semi-labile DOC 1/250days
       DON_reminR  = (c1/160.0_r8) * dps,          & ! rate for semi-labile DON 1/160days
       DOFe_reminR = (c1/160.0_r8) * dps,          & ! rate for semi-labile DOFe 1/160days
       DOP_reminR  = (c1/160.0_r8) * dps,          & ! rate for semi-labile DOP 1/160days  
       DONr_reminR = (c1/(365.0_r8*2.5_r8)) * dps, & ! timescale for refrac DON 1/2.5yrs
       DOPr_reminR = (c1/(365.0_r8*2.5_r8)) * dps, & ! timescale for refrac DOP 1/2.5yrs
       DONrefract = 0.08_r8,                       & ! fraction of DON to refractory pool
       DOPrefract = 0.03_r8                          ! fraction of DOP to refractory pool

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE ecosys_parms_init

    USE io_types, ONLY: stdout, nml_in, nml_filename
    USE broadcast, ONLY : broadcast_scalar, broadcast_array

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: subname = 'ecosys_parms:ecosys_parms_init'

    INTEGER(KIND=int_kind) :: auto_ind

    LOGICAL(KIND=log_kind) :: &
         lnml_found             ! Was ecosys_parms_nml found ?

    NAMELIST /ecosys_parms_nml/ &
         parm_Fe_bioavail, &
         parm_o2_min, &
         parm_o2_min_delta, &
         parm_kappa_nitrif, &
         parm_nitrif_par_lim, &
         parm_z_mort_0, &
         parm_z_mort2_0, &
         parm_labile_ratio, &
         parm_POMbury, &
         parm_BSIbury, &
         parm_fe_scavenge_rate0, &
         parm_f_prod_sp_CaCO3, &
         parm_POC_diss, &
         parm_SiO2_diss, &
         parm_CaCO3_diss, &
         parm_scalelen_z, &
         parm_scalelen_vals, &
         autotrophs

    !---------------------------------------------------------------------------
    !   default namelist settings
    !---------------------------------------------------------------------------

    parm_Fe_bioavail       = 1.0_r8
    parm_o2_min            = 4.0_r8
    parm_o2_min_delta      = 2.0_r8
    parm_kappa_nitrif      = 0.06_r8 * dps  ! (= 1/( days))
    parm_nitrif_par_lim    = 1.0_r8
    parm_z_mort_0          = 0.1_r8 * dps
    parm_z_mort2_0         = 0.4_r8 * dps
    parm_labile_ratio      = 0.85_r8
    parm_POMbury           = 1.4_r8         ! x1 default
    parm_BSIbury           = 0.65_r8        ! x1 default
    parm_fe_scavenge_rate0 = 3.0_r8         ! x1 default
    parm_f_prod_sp_CaCO3   = 0.055_r8       ! x1 default
    parm_POC_diss          = 88.0e2_r8
    parm_SiO2_diss         = 250.0e2_r8
    parm_CaCO3_diss        = 150.0e2_r8

    parm_scalelen_z    = (/ 130.0e2_r8, 290.0e2_r8, 670.0e2_r8, 1700.0e2_r8 /) ! x1 default
    parm_scalelen_vals = (/     1.0_r8,     3.0_r8,     5.0_r8,      9.0_r8 /) ! x1 default

    auto_ind = sp_ind
    autotrophs(auto_ind)%sname         = 'sp'
    autotrophs(auto_ind)%lname         = 'Small Phyto'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .true.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = auto_ind
    autotrophs(auto_ind)%kFe           = 0.04e-3_r8
    autotrophs(auto_ind)%kPO4          = 0.01_r8
    autotrophs(auto_ind)%kDOP          = 0.26_r8
    autotrophs(auto_ind)%kNO3          = 0.1_r8
    autotrophs(auto_ind)%kNH4          = 0.01_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_r8
    autotrophs(auto_ind)%Qp            = 0.00855_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_r8
    autotrophs(auto_ind)%alphaPI       = 0.6_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_r8
    autotrophs(auto_ind)%temp_thres    = -10.0_r8
    autotrophs(auto_ind)%mort          = 0.12_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.01_r8
    autotrophs(auto_ind)%z_umax_0      = 3.3_r8 * dps ! x1 default
    autotrophs(auto_ind)%z_grz         = 1.05_r8              
    autotrophs(auto_ind)%graze_zoo     = 0.3_r8
    autotrophs(auto_ind)%graze_poc     = 0.0_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.15_r8

    auto_ind = diat_ind
    autotrophs(auto_ind)%sname         = 'diat'
    autotrophs(auto_ind)%lname         = 'Diatom'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = auto_ind
    autotrophs(auto_ind)%kFe           = 0.06e-3_r8
    autotrophs(auto_ind)%kPO4          = 0.05_r8
    autotrophs(auto_ind)%kDOP          = 0.9_r8
    autotrophs(auto_ind)%kNO3          = 0.5_r8
    autotrophs(auto_ind)%kNH4          = 0.05_r8
    autotrophs(auto_ind)%kSiO3         = 0.8_r8
    autotrophs(auto_ind)%Qp            = 0.00855_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_r8
    autotrophs(auto_ind)%alphaPI       = 0.465_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 4.0_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_r8
    autotrophs(auto_ind)%temp_thres    = -10.0_r8
    autotrophs(auto_ind)%mort          = 0.12_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.02_r8
    autotrophs(auto_ind)%z_umax_0      = 3.08_r8 * dps ! x1 default
    autotrophs(auto_ind)%z_grz         = 1.0_r8              
    autotrophs(auto_ind)%graze_zoo     = 0.3_r8
    autotrophs(auto_ind)%graze_poc     = 0.42_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.2_r8

    auto_ind = diaz_ind
    autotrophs(auto_ind)%sname         = 'diaz'
    autotrophs(auto_ind)%lname         = 'Diazotroph'
    autotrophs(auto_ind)%Nfixer        = .true.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = auto_ind
    autotrophs(auto_ind)%kFe           = 0.04e-3_r8
    autotrophs(auto_ind)%kPO4          = 0.02_r8
    autotrophs(auto_ind)%kDOP          = 0.09_r8
    autotrophs(auto_ind)%kNO3          = 1.0_r8
    autotrophs(auto_ind)%kNH4          = 0.15_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_r8
    autotrophs(auto_ind)%Qp            = 0.002735_r8
    autotrophs(auto_ind)%gQfe_0        = 60.0e-6_r8
    autotrophs(auto_ind)%gQfe_min      = 12.0e-6_r8
    autotrophs(auto_ind)%alphaPI       = 0.4_r8 * dps
    autotrophs(auto_ind)%PCref         = 0.7_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_r8
    autotrophs(auto_ind)%loss_thres    = 0.022_r8
    autotrophs(auto_ind)%loss_thres2   = 0.001_r8
    autotrophs(auto_ind)%temp_thres    = 14.0_r8
    autotrophs(auto_ind)%mort          = 0.15_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.0_r8
    autotrophs(auto_ind)%agg_rate_max  = 0.0_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.0_r8
    autotrophs(auto_ind)%z_umax_0      = 0.6_r8 * dps
    autotrophs(auto_ind)%z_grz         = 1.2_r8              
    autotrophs(auto_ind)%graze_zoo     = 0.3_r8
    autotrophs(auto_ind)%graze_poc     = 0.05_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.15_r8

    !---------------------------------------------------------------------------
    !   read in namelist
    !---------------------------------------------------------------------------

    IF (my_task == master_task) THEN
       lnml_found = .FALSE.
       OPEN(UNIT=nml_in, FILE=nml_filename, STATUS='OLD')
10     CONTINUE
       READ(UNIT=nml_in, NML=ecosys_parms_nml, ERR=10, END=20)
       CLOSE(UNIT=nml_in)
       lnml_found = .TRUE.
20     CONTINUE
    END IF

    CALL broadcast_scalar(lnml_found, master_task)
    IF (.NOT. lnml_found) THEN
       CALL document(subname, 'ecosys_parms_nml not found')
       CALL exit_POP(sigAbort, 'ERROR : stopping in ' // subname)
    END IF

    !---------------------------------------------------------------------------
    !   broadcast all namelist variables
    !---------------------------------------------------------------------------

    CALL broadcast_scalar(parm_Fe_bioavail, master_task)
    CALL broadcast_scalar(parm_o2_min, master_task)
    CALL broadcast_scalar(parm_o2_min_delta, master_task)
    CALL broadcast_scalar(parm_kappa_nitrif, master_task)
    CALL broadcast_scalar(parm_nitrif_par_lim, master_task)
    CALL broadcast_scalar(parm_z_mort_0, master_task)
    CALL broadcast_scalar(parm_z_mort2_0, master_task)
    CALL broadcast_scalar(parm_labile_ratio, master_task)
    CALL broadcast_scalar(parm_POMbury, master_task)
    CALL broadcast_scalar(parm_BSIbury, master_task)
    CALL broadcast_scalar(parm_fe_scavenge_rate0, master_task)
    CALL broadcast_scalar(parm_f_prod_sp_CaCO3, master_task)
    CALL broadcast_scalar(parm_POC_diss, master_task)
    CALL broadcast_scalar(parm_SiO2_diss, master_task)
    CALL broadcast_scalar(parm_CaCO3_diss, master_task)

    CALL broadcast_array(parm_scalelen_z, master_task)
    CALL broadcast_array(parm_scalelen_vals, master_task)

    DO auto_ind = 1, autotroph_cnt
       CALL broadcast_scalar(autotrophs(auto_ind)%sname, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%lname, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%Nfixer, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%imp_calcifier, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%exp_calcifier, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%grazee_ind, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%kFe, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%kPO4, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%kDOP, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%kNO3, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%kNH4, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%kSiO3, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%Qp, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%gQfe_0, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%gQfe_min, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%alphaPI, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%PCref, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%thetaN_max, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%loss_thres, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%loss_thres2, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%temp_thres, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%mort, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%mort2, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%agg_rate_max, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%agg_rate_min, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%z_umax_0, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%z_grz, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%graze_zoo, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%graze_poc, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%graze_doc, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%loss_poc, master_task)
       CALL broadcast_scalar(autotrophs(auto_ind)%f_zoo_detr, master_task)
    END DO

    !---------------------------------------------------------------------------
    !   echo all namelist variables to stdout
    !---------------------------------------------------------------------------

    IF (my_task == master_task) THEN
       WRITE (stdout,*) '----------------------------------------'
       WRITE (stdout,*) '----- ecosys_parms namelist values -----'
       WRITE (stdout,*) 'parm_Fe_bioavail       = ', parm_Fe_bioavail
       WRITE (stdout,*) 'parm_o2_min            = ', parm_o2_min
       WRITE (stdout,*) 'parm_o2_min_delta      = ', parm_o2_min_delta
       WRITE (stdout,*) 'parm_kappa_nitrif      = ', parm_kappa_nitrif
       WRITE (stdout,*) 'parm_nitrif_par_lim    = ', parm_nitrif_par_lim
       WRITE (stdout,*) 'parm_z_mort_0          = ', parm_z_mort_0
       WRITE (stdout,*) 'parm_z_mort2_0         = ', parm_z_mort2_0
       WRITE (stdout,*) 'parm_labile_ratio      = ', parm_labile_ratio
       WRITE (stdout,*) 'parm_POMbury           = ', parm_POMbury
       WRITE (stdout,*) 'parm_BSIbury           = ', parm_BSIbury
       WRITE (stdout,*) 'parm_fe_scavenge_rate0 = ', parm_fe_scavenge_rate0
       WRITE (stdout,*) 'parm_f_prod_sp_CaCO3   = ', parm_f_prod_sp_CaCO3
       WRITE (stdout,*) 'parm_POC_diss          = ', parm_POC_diss
       WRITE (stdout,*) 'parm_SiO2_diss         = ', parm_SiO2_diss
       WRITE (stdout,*) 'parm_CaCO3_diss        = ', parm_CaCO3_diss
       WRITE (stdout,*) 'parm_scalelen_z        = ', parm_scalelen_z
       WRITE (stdout,*) 'parm_scalelen_vals     = ', parm_scalelen_vals
       DO auto_ind = 1, autotroph_cnt
          WRITE (stdout,*) 'lname(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%lname
          WRITE (stdout,*) 'Nfixer(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Nfixer
          WRITE (stdout,*) 'imp_calcifier(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%imp_calcifier
          WRITE (stdout,*) 'exp_calcifier(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%exp_calcifier
          WRITE (stdout,*) 'grazee_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%grazee_ind
          WRITE (stdout,*) 'kFe(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kFe
          WRITE (stdout,*) 'kPO4(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kPO4
          WRITE (stdout,*) 'kDOP(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kDOP
          WRITE (stdout,*) 'kNO3(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kNO3
          WRITE (stdout,*) 'kNH4(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kNH4
          WRITE (stdout,*) 'kSiO3(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kSiO3
          WRITE (stdout,*) 'Qp(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Qp
          WRITE (stdout,*) 'gQfe_0(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%gQfe_0
          WRITE (stdout,*) 'gQfe_min(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%gQfe_min
          WRITE (stdout,*) 'alphaPI(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%alphaPI
          WRITE (stdout,*) 'PCref(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%PCref
          WRITE (stdout,*) 'thetaN_max(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%thetaN_max
          WRITE (stdout,*) 'loss_thres(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%loss_thres
          WRITE (stdout,*) 'loss_thres2(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%loss_thres2
          WRITE (stdout,*) 'temp_thres(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%temp_thres
          WRITE (stdout,*) 'mort(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%mort
          WRITE (stdout,*) 'mort2(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%mort2
          WRITE (stdout,*) 'agg_rate_max(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%agg_rate_max
          WRITE (stdout,*) 'agg_rate_min(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%agg_rate_min
          WRITE (stdout,*) 'z_umax_0(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%z_umax_0
          WRITE (stdout,*) 'z_grz(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%z_grz
          WRITE (stdout,*) 'graze_zoo(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%graze_zoo
          WRITE (stdout,*) 'graze_poc(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%graze_poc
          WRITE (stdout,*) 'graze_doc(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%graze_doc
          WRITE (stdout,*) 'loss_poc(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%loss_poc
          WRITE (stdout,*) 'f_zoo_detr(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%f_zoo_detr
       END DO
       WRITE (stdout,*) '----------------------------------------'
    END IF

  END SUBROUTINE ecosys_parms_init

  !*****************************************************************************

END MODULE ecosys_parms
