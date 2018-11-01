#ifdef CLM_PFLOTRAN

module clm_pflotran_interface_data
!
! NOTES for convenience:
!        (1) '*_pfp': mpi vecs for PF variables; '_clmp': mpi vecs for CLM variables;
!            '*_pfs': seq vecs for PF variables; '_clms': seq. vecs for CLM variables;
!        (2) '*_': 3D (XYZ) subsurface domain's variables;
!                  with '_x/y/z_' used for different directions of 3D domain.
!        (3) '*_subsurf_': 2D (XY) surface of 3D domain's varialbes;
!                         for bottom, uses '_subbase_', but essentially supposing same shape/area as surface.
!            '*_srf_': for variables with 2D-grids at ground surface, which may or may not same as '_subsurf_'. NOTE that it's not supported now.
! Revised by Fengming Yuan, CCSI-ORNL @May-2015

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
  use petscsys
  use petscvec

  implicit none

  private

  type, public :: clm_pflotran_idata_type


  !------------------- A few global constants --------------------------------------------
  
  ! numbers of CLM soil layers, grids that are mapped to/from PFLOTRAN (global constants, not local copy)
  PetscInt :: nzclm_mapped
  PetscInt :: nxclm_mapped
  PetscInt :: nyclm_mapped
  PetscReal :: x0clm_global
  PetscReal :: y0clm_global
  PetscReal :: z0clm_global
  PetscReal, pointer :: dxclm_global(:)              ! this is NOT the 3-D vec 'dxsoil' defined below, rather it's the universal x-direction interval (OR, longitudal degree interval from CLM land surf grids) for all gridcells
  PetscReal, pointer :: dyclm_global(:)              ! this is NOT the 3-D vec 'dysoil' defined below, rather it's the universal y-direction interval (OR, longitudal degree interval from CLM land surf grids)
  PetscReal, pointer :: dzclm_global(:)              ! this is NOT the 3-D vec 'dzsoil' defined below, rather it's the universal soil layer thickness (unit: m) for all gridcells

  ! --- Decompose domain in 3-D (only work with structured PF grid currently) ------------
  
  ! processors no.
  PetscInt :: npx, npy, npz
  ! domain nodes no. for each processors
  PetscInt, pointer :: clm_lx(:)   ! array size is 'npx'
  PetscInt, pointer :: clm_ly(:)   ! array size is 'npy'
  PetscInt, pointer :: clm_lz(:)   ! array size is 'npz'

  !--------------------  Mesh property ---------------------------------------------------
  
  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_sub    ! num of local clm cells
  PetscInt :: ngclm_sub    ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub     ! num of local pflotran cells
  PetscInt :: ngpf_sub     ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the top/bottom cells of the 3D subsurface domain
  PetscInt :: nlclm_2dtop  ! num of local clm cells
  PetscInt :: ngclm_2dtop  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dtop   ! num of local pflotran cells
  PetscInt :: ngpf_2dtop   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  PetscInt :: nlclm_2dbot  ! num of local clm cells
  PetscInt :: ngclm_2dbot  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dbot   ! num of local pflotran cells
  PetscInt :: ngpf_2dbot   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlclm_srf    ! num of local clm cells
  PetscInt :: ngclm_srf    ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf     ! num of local pflotran cells
  PetscInt :: ngpf_srf     ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! sub-surf/sub-base area (in 2D): the toppest/lowest cells of subsurface domain for BCs
  ! At this moment, assumes both surf/base cells are exactly SAME in area.
  Vec :: area_subsurf_clmp   ! mpi vec
  Vec :: area_subsurf_pfs    ! seq vec
  Vec :: area_subsurf_pfp    ! mpi vec
  Vec :: area_subsurf_clms   ! seq vec

  ! Area of top face of PF cells (in 3D) (note: a PF cell has faces (facets) of top/bottom/east/west/south/north)
  Vec :: area_top_face_clmp  ! mpi vec
  Vec :: area_top_face_pfs   ! seq vec
  Vec :: area_top_face_pfp   ! mpi vec
  Vec :: area_top_face_clms  ! seq vec

  ! z-axis (in 3D), soil depth of center of a soil cell in unit of meters
  Vec :: zsoil_clmp          ! mpi vec
  Vec :: zsoil_pfs           ! seq vec
  Vec :: zsoil_pfp           ! mpi vec
  Vec :: zsoil_clms          ! seq vec

  ! x/y-axis (in 3D), grid center of a soil cell in unit of meters
  Vec :: xsoil_clmp          ! mpi vec
  Vec :: xsoil_pfs           ! seq vec
  Vec :: ysoil_clmp          ! mpi vec
  Vec :: ysoil_pfs           ! seq vec

  ! soil cell inter-nodes coordinates ('vertex' called in PF mesh; 'interface level' called in CLM soil layers)
  Vec :: zisoil_clmp          ! mpi vec
  Vec :: zisoil_pfs           ! seq vec

  ! length/width/thickness of soil cells (in 3D) in unit of meters
  Vec :: dxsoil_clmp          ! mpi vec
  Vec :: dxsoil_pfs           ! seq vec
  Vec :: dysoil_clmp          ! mpi vec
  Vec :: dysoil_pfs           ! seq vec
  Vec :: dzsoil_clmp          ! mpi vec
  Vec :: dzsoil_pfs           ! seq vec
  ! a NOTE here: Given a 3D-cell's 'area_top_face' and 'zsoi' known, 
  !              it's possible to calculate its volume (may be useful ?)

   ! cell IDs (in 3D) (for tesing meshes)
  Vec :: cellid_clmp          ! mpi vec
  Vec :: cellid_pfs           ! seq vec
  Vec :: cellid_pfp           ! mpi vec
  Vec :: cellid_clms          ! seq vec
   ! top layer cell IDs (in 2D) (for tesing meshes)
  Vec :: cellid_2dtop_clmp    ! mpi vec
  Vec :: cellid_2dtop_pfs     ! seq vec
  Vec :: cellid_2dtop_pfp     ! mpi vec
  Vec :: cellid_2dtop_clms    ! seq vec

  !-------------- TH properties ----------------------------------------------------------
  
  PetscBool :: head_based
  PetscReal :: pressure_reference

  ! CLM's hydraulic properties
  Vec :: hksat_x_clmp
  Vec :: hksat_y_clmp
  Vec :: hksat_z_clmp
  Vec :: watsat_clmp
  Vec :: watfc_clmp
  Vec :: bulkdensity_dry_clmp
  Vec :: effporosity_clmp

  Vec :: hksat_x_pfs
  Vec :: hksat_y_pfs
  Vec :: hksat_z_pfs
  Vec :: watsat_pfs
  Vec :: watfc_pfs
  Vec :: bulkdensity_dry_pfs
  Vec :: effporosity_pfs

  ! clapp-Horburger's function parameters - needed in BGC somehow
  Vec :: sucsat_clmp   
  Vec :: bsw_clmp
  Vec :: sucsat_pfs
  Vec :: bsw_pfs

  ! PF's hydraulic properties
  Vec :: sr_pcwmax_pfp
  Vec :: pcwmax_pfp
  Vec :: effporosity_pfp
  Vec :: sr_pcwmax_clms
  Vec :: pcwmax_clms
  Vec :: effporosity_clms

  ! CLM's thermal properties
  Vec :: tkwet_clmp     ! unit: W/m/K
  Vec :: tkdry_clmp
  Vec :: tkfrz_clmp
  Vec :: hcvsol_clmp    ! unit: J/m^3-K

  Vec :: tkwet_pfs
  Vec :: tkdry_pfs
  Vec :: tkfrz_pfs
  Vec :: hcvsol_pfs

  ! TH state vecs from CLM (mpi) to PF (seq)
  Vec :: press_ref_clmp                 ! reference pressure head (Pa)
  Vec :: press_ref_pfs

  Vec :: press_clmp                     ! water pressure head (Pa)
  Vec :: soilpsi_clmp                   ! soil matric potential (Pa)
  Vec :: soillsat_clmp                  ! soil liq. water saturation (0 - 1)
  Vec :: soilisat_clmp                  ! soil ice water saturation (0 - 1)
  Vec :: soilliq_clmp                   ! soil liq. water mass (kg/m3 bulk soil)
  Vec :: soilice_clmp                   ! soil ice water mass (kg/m3 bulk soil)
  Vec :: soilt_clmp                     ! soil temperature (degC)
  Vec :: press_pfs
  Vec :: soilpsi_pfs
  Vec :: soillsat_pfs
  Vec :: soilisat_pfs
  Vec :: soilliq_pfs
  Vec :: soilice_pfs
  Vec :: soilt_pfs

  ! TH state vecs from PF (mpi) to CLM (seq)

  Vec :: press_pfp                     ! water pressure head (Pa)
  Vec :: soilpsi_pfp                   ! soil matric potential (Pa)
  Vec :: soillsat_pfp                  ! soil liq. water saturation (0 - 1)
  Vec :: soilisat_pfp                  ! soil ice water saturation (0 - 1)
  Vec :: soilliq_pfp                   ! soil liq. water mass (kg/m3 bulk soil)
  Vec :: soilice_pfp                   ! soil ice water mass (kg/m3 bulk soil)
  Vec :: soilt_pfp                     ! soil temperature (degC)
  Vec :: press_clms
  Vec :: soilpsi_clms
  Vec :: soillsat_clms
  Vec :: soilisat_clms
  Vec :: soilliq_clms
  Vec :: soilice_clms
  Vec :: soilt_clms
 
  !------------------------------- Soil BGC ----------------------------------------------

  !
  ! ----- BGC constants
  !

  ! the following constants are for consistently coverting mass and moles btw CLM-CN and PF bgc
  ! so that mass balance error would not be caused by these
  PetscReal :: N_molecular_weight
  PetscReal :: C_molecular_weight

  ! Soil BGC decomposing pools
  PetscInt :: ndecomp_pools
  logical, pointer :: floating_cn_ratio(:)           ! TRUE => pool has variable C:N ratio
  PetscInt :: ndecomp_elements                       ! no. of elements considered: C, N
  PetscReal, pointer:: decomp_element_ratios(:,:)    ! ratios of elements in decomposing pools (unit: moles)

  !NOTES: The following is what PF bgc right now using for CLM-PFLOTRAN coupling
  ! if need adding or modifying, it's possible and update BOTH here and subroutine 'pflotranModelGetRTspecies'
  ! (Of course, it must be modifying the PF input card and get those variables and relevant reactions in RT).

  ! RT bgc species 'idof' and 'name'
  PetscInt, pointer:: ispec_decomp_c(:)              ! name: pool_name, OR, pool_name // "C"
  PetscInt, pointer:: ispec_decomp_n(:)              ! name: "", OR, pool_name // "N"
  PetscInt, pointer:: ispec_decomp_hr(:)             ! name: pool_name // "CHR"
  PetscInt, pointer:: ispec_decomp_nmin(:)           ! name: pool_name // "NMIN"
  PetscInt, pointer:: ispec_decomp_nimp(:)           ! name: pool_name // "NIMM"
  PetscInt, pointer:: ispec_decomp_nimm(:)           ! name: pool_name // "NIMP"
  character(len=32), pointer :: decomp_pool_name(:)  ! name of pools
  PetscReal, pointer:: ck_decomp_c(:)                ! K: first-order decomposition rate constant in 1/sec
  PetscReal, pointer:: adfactor_ck_c(:)              ! scalar to adjust K based on decomp pools - currently used for speeding-up decomposition
  PetscReal, pointer:: fr_decomp_c(:,:)              ! fractions of downstream pools (receivers - pools excluding donor but adding CO2, note (k,k) is that fraction of CO2 respired)

  PetscInt:: ispec_hrimm
  character(len=32):: name_hrim   = "HRimm"          ! this is for total HR

  PetscInt :: ispec_nmin, ispec_nimp, ispec_nimm
  character(len=32):: name_nmin  = "Nmin"            ! this is for total Nmin
  character(len=32):: name_nimp  = "Nimp"            ! this is for total Nimmp
  character(len=32):: name_nimm  = "Nimm"            ! this is for total Nimm

  PetscInt:: ispec_nh4, ispec_no3, ispec_nh4s, ispec_no3s, ispec_nh4sorb
  character(len=32):: name_nh4     = "NH4+"
  character(len=32):: name_no3     = "NO3-"
  character(len=32):: name_nh4s    = "Ammonium"
  character(len=32):: name_no3s    = "Nitrate"
  character(len=32):: name_nh4sorb = "NH4sorb"

  PetscInt :: ispec_plantndemand, ispec_plantnh4uptake, ispec_plantno3uptake
  character(len=32):: name_plantndemand   = "Plantndemand"
  character(len=32):: name_plantnh4uptake = "Plantnh4uptake"
  character(len=32):: name_plantno3uptake = "Plantno3uptake"

  PetscInt :: ispec_ngasmin, ispec_ngasnitr, ispec_ngasdeni
  character(len=32):: name_ngasmin = "NGASmin"
  character(len=32):: name_ngasnitr= "NGASnitr"
  character(len=32):: name_ngasdeni= "NGASdeni"

  PetscInt :: ispec_co2aq, ispec_n2aq, ispec_n2oaq
  character(len=32):: name_co2aq   = "CO2(aq)"
  character(len=32):: name_n2aq    = "N2(aq)"
  character(len=32):: name_n2oaq   = "N2O(aq)"

  PetscInt :: ispec_co2, ispec_n2, ispec_n2o   ! this is for gases as a holder (immobile)
  character(len=32):: name_co2 = "CO2imm"
  character(len=32):: name_n2o = "N2Oimm"
  character(len=32):: name_n2  = "N2imm"

  !
  ! -----BGC vecs from CLM (mpi, ghosted) to PF (seq, local)
  !
  ! the following can be directly used to drive BGC
  Vec :: t_scalar_clmp                  ! temperature response function value from CLM for decomposition reations
  Vec :: w_scalar_clmp                  ! soil moisture response function value from CLM for decomposition reations
  Vec :: o_scalar_clmp                  ! soil anoxic response function value from CLM for decomposition and/or CH4/NOx processes
  Vec :: depth_scalar_clmp              ! soil depth response function value from CLM for decomposition reations
  Vec :: t_scalar_pfs
  Vec :: w_scalar_pfs
  Vec :: o_scalar_pfs
  Vec :: depth_scalar_pfs
  
  ! A NOTE here: the folllowing decomposing pool vec (1D) is ordered by 'cell' first, then 'species'

  ! initial ground/soil C/N pools from CLM (mpi) to PF (seq)
  Vec :: decomp_cpools_vr_clmp     ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_clmp     ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_clmp          ! (moleN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_clmp          ! (moleN/m3) vertically-resolved soil mineral NH4
  Vec :: smin_nh4sorb_vr_clmp      ! (moleN/m3) vertically-resolved soil absorbed mineral NH4
  Vec :: decomp_cpools_vr_pfs      ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_pfs      ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_pfs           ! (moleN/m3) vertically-resolved soil mineral NO3
  Vec :: smin_nh4_vr_pfs           ! (moleN/m3) vertically-resolved soil mineral NH4
  Vec :: smin_nh4sorb_vr_pfs       ! (moleN/m3) vertically-resolved soil absorbed mineral NH4

  ! time-varying ground/soil C/N rates from CLM (mpi) to PF (seq) (Units: moleC(N)/m3/s)
  Vec :: kscalar_decomp_c_clmp     ! (unitless) a site scalar (c,j) to adjust SOM decomposition rate constants (default: 1.0)
  Vec :: rate_decomp_c_clmp
  Vec :: rate_decomp_n_clmp
  Vec :: rate_smin_no3_clmp
  Vec :: rate_smin_nh4_clmp
  Vec :: rate_plantndemand_clmp
  Vec :: kscalar_decomp_c_pfs
  Vec :: rate_decomp_c_pfs
  Vec :: rate_decomp_n_pfs
  Vec :: rate_smin_no3_pfs
  Vec :: rate_smin_nh4_pfs
  Vec :: rate_plantndemand_pfs

  !
  ! -----BGC vecs from PF (mpi, ghosted) to CLM (seq, local)
  !
  ! BGC state variables
  Vec :: decomp_cpools_vr_pfp
  Vec :: decomp_npools_vr_pfp
  Vec :: smin_no3_vr_pfp
  Vec :: smin_nh4_vr_pfp                ! (moleN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_vr_pfp            ! (moleN/m3) vertically-resolved absorbed NH4-N
  !
  Vec :: decomp_cpools_vr_clms          ! (moleC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
  Vec :: decomp_npools_vr_clms          ! (moleN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
  Vec :: smin_no3_vr_clms               ! (moleN/m3) vertically-resolved total soil mineral NO3
  Vec :: smin_nh4_vr_clms               ! (moleN/m3) vertically-resolved total soil mineral NH4 (incl. absorbed)
  Vec :: smin_nh4sorb_vr_clms           ! (moleN/m3) vertically-resolved absorbed NH4-N
  !
  ! 'accextrn' is accumulative N extract by plant roots in 'PFLOTRAN' within a CLM timestep
  Vec :: accextrnh4_vr_pfp                ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrnh4_vr_clms               ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrno3_vr_pfp                ! (moleN/m3) vertically-resolved root extraction N
  Vec :: accextrno3_vr_clms               ! (moleN/m3) vertically-resolved root extraction N

  ! gases in water (aqueous solution of gases)
  ! gases species is accumulative in 'PFLOTRAN', so needs to calculate their fluxes in the CLM-PF interface and reset back to PFLOTRAN
  Vec :: gco2_vr_pfp                   ! (moleC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_vr_clms                  ! (moleC/m3) vertically-resolved soil CO2 C
  Vec :: gco2_vr_clmp                  ! (moleC/m3) vertically-resolved soil CO2 C, after gas emission
  Vec :: gco2_vr_pfs                   ! (moleC/m3) vertically-resolved soil CO2 C, after gas emission

  Vec :: gn2_vr_pfp                    ! (moleN/m3) vertically-resolved N2-N
  Vec :: gn2_vr_clms                   ! (moleN/m3) vertically-resolved N2-N
  Vec :: gn2_vr_clmp                   ! (moleN/m3) vertically-resolved N2-N, after gas emission
  Vec :: gn2_vr_pfs                    ! (moleN/m3) vertically-resolved N2-N, after gas emission

  Vec :: gn2o_vr_pfp                   ! (moleN/m3) vertically-resolved N2O-N
  Vec :: gn2o_vr_clms                  ! (moleN/m3) vertically-resolved N2O-N
  Vec :: gn2o_vr_clmp                  ! (moleN/m3) vertically-resolved N2O-N, after gas emission
  Vec :: gn2o_vr_pfs                   ! (moleN/m3) vertically-resolved N2O-N, after gas emission

  ! some tracking variables from PFLOTRAN bgc to obtain reaction flux rates which needed by CLM
  Vec :: acchr_vr_pfp                  ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from individual decomposition
  Vec :: acchr_vr_clms                 ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from individual decomposition
  Vec :: acctothr_vr_pfp               ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from all decomposition
  Vec :: acctothr_vr_clms              ! (moleC/m3/timestep) vertically-resolved heterotrophic resp. C from all decomposition

  Vec :: accnmin_vr_pfp                ! (moleN/m3/timestep) vertically-resolved N mineralization
  Vec :: accnmin_vr_clms               ! (moleN/m3/timestep) vertically-resolved N mineralization
  Vec :: acctotnmin_vr_pfp             ! (moleN/m3/timestep) vertically-resolved total N mineralization
  Vec :: acctotnmin_vr_clms            ! (moleN/m3/timestep) vertically-resolved total N mineralization

  Vec :: accnimmp_vr_pfp               ! (moleN/m3/timestep) vertically-resolved potential N immoblization
  Vec :: accnimmp_vr_clms              ! (moleN/m3/timestep) vertically-resolved potential N immoblization
  Vec :: acctotnimmp_vr_pfp            ! (moleN/m3/timestep) vertically-resolved potential total N immoblization
  Vec :: acctotnimmp_vr_clms           ! (moleN/m3/timestep) vertically-resolved potential total N immoblization

  Vec :: accnimm_vr_pfp                ! (moleN/m3/timestep) vertically-resolved N immoblization
  Vec :: accnimm_vr_clms               ! (moleN/m3/timestep) vertically-resolved N immoblization
  Vec :: acctotnimm_vr_pfp             ! (moleN/m3/timestep) vertically-resolved total N immoblization
  Vec :: acctotnimm_vr_clms            ! (moleN/m3/timestep) vertically-resolved total N immoblization

  Vec :: accngasmin_vr_pfp              ! (moleN/m3/timestep) vertically-resolved N2O-N from mineralization
  Vec :: accngasmin_vr_clms             ! (moleN/m3/timestep) vertically-resolved N2O-N from mineralization

  Vec :: accngasnitr_vr_pfp             ! (moleN/m3/timestep) vertically-resolved N2O-N from nitrification
  Vec :: accngasnitr_vr_clms            ! (moleN/m3/timestep) vertically-resolved N2O-N from nitrification

  Vec :: accngasdeni_vr_pfp             ! (moleN/m3/timestep) vertically-resolved N2-N from denitrification
  Vec :: accngasdeni_vr_clms            ! (moleN/m3/timestep) vertically-resolved N2-N from denitrification

  ! actual aqeuous N mass flow rate(moleN/m2/sec) at the top (runoff)/bottom (leaching) of 3-D subsurface domain
  ! (+ in, - out)
  Vec :: f_nh4_subsurf_pfp    ! mpi vec
  Vec :: f_nh4_subsurf_clms   ! seq vec
  Vec :: f_nh4_subbase_pfp    ! mpi vec
  Vec :: f_nh4_subbase_clms   ! seq vec
  Vec :: f_no3_subsurf_pfp    ! mpi vec
  Vec :: f_no3_subsurf_clms   ! seq vec
  Vec :: f_no3_subbase_pfp    ! mpi vec
  Vec :: f_no3_subbase_clms   ! seq vec


  !------------------------------- Soil Thermal-Hydrology ----------------------------------------------

  !
  ! -----TH vecs from CLM (mpi, ghosted) to PF (seq, local)
  !

  ! Sink/Source of water (with thermal) for PFLOTRAN's 3D subsurface domain
  Vec :: qflow_clmp   ! mpi vec (H2O): kgH2O/m3/sec
  Vec :: qflow_pfs    ! seq vec
  Vec :: qflowt_clmp  ! mpi vec (H2O)
  Vec :: qflowt_pfs   ! seq vec

  ! Sink/Source of (non-mass) heat flow rate for 3-D subsurface domain
  Vec :: eflow_clmp   ! mpi vec: all non-mass forms of energy src/sink rate (MJ/m3/sec)
  Vec :: eflow_pfs    ! seq vec

  ! BC: water pressure (Pa) on the 2D top/bottom interface of 3-D subsurface domain as boundary conditions from CLM to PF
  Vec :: press_subsurf_clmp    ! mpi vec
  Vec :: press_subbase_clmp    ! mpi vec
  Vec :: press_subsurf_pfs     ! seq vec
  Vec :: press_subbase_pfs     ! seq vec

  Vec :: press_maxponding_clmp ! mpi vec
  Vec :: press_maxponding_pfs  ! seq vec

  ! BC-h: water infiltration/recharge(drainage) (mH2O/sec) on the 2D top/bottom interface of 3-D subsurface domain as boundary conditions from CLM to PF
  Vec :: qfluxw_subsurf_clmp    ! mpi vec
  Vec :: qfluxw_subbase_clmp    ! mpi vec
  Vec :: qfluxw_subsurf_pfs     ! seq vec
  Vec :: qfluxw_subbase_pfs     ! seq vec

  ! BC-h: soil evaporation (mH2O/sec) on the 2D top interface of 3-D subsurface domain as boundary conditions from CLM to PF
  Vec :: qfluxev_subsurf_clmp   ! mpi vec
  Vec :: qfluxev_subsurf_pfs    ! seq vec

  ! BC-t: temperature/eflux at the subsurface interface (2D TOP)
  Vec :: gtemp_subsurf_clmp   ! mpi vec: (1) for specifying temperature (thermal state, like enthalpy), or (2) for heat conductance at BC
  Vec :: gtemp_subsurf_pfs    ! seq vec
  Vec :: eflux_subsurf_clmp   ! mpi vec: all forms of energy (MJ/m^2-sec)
  Vec :: eflux_subsurf_pfs    ! seq vec

  Vec :: efluxr_subsurf_clmp  ! mpi vec: radiation form of energy (MJ/m^2-sec)
  Vec :: efluxr_subsurf_pfs   ! seq vec

  Vec :: efluxl_subsurf_clmp  ! mpi vec: latent heat form of energy (MJ/m^2-sec)
  Vec :: efluxl_subsurf_pfs   ! seq vec
  ! (a note here: if specifying 'gtemp_subsurf' above, sensible heat flux at interface may not be needed)

  ! BC-t: temperature/eflux at the subsurface interface (2D BOTTOM)
  Vec :: eflux_subbase_clmp  ! mpi vec
  Vec :: eflux_subbase_pfs   ! seq vec
  Vec :: gtemp_subbase_clmp  ! mpi vec
  Vec :: gtemp_subbase_pfs   ! seq vec

  !
  ! -----TH vecs from PF (mpi, ghosted) to CLM (seq, local)
  !

  ! actual mass water flow rate (kgH2O/sec) through the top/bottom BC (2-D) of 3-D subsurface domain
  ! (+ in, - out)
  Vec :: qinfl_subsurf_pfp    ! mpi vec: actual infiltration (+)
  Vec :: qinfl_subsurf_clms   ! seq vec
  Vec :: qsurf_subsurf_pfp    ! mpi vec: actual overland flow - potential-actual infiltration or water upwarding (-)
  Vec :: qsurf_subsurf_clms   ! seq vec
  Vec :: qflux_subbase_pfp    ! mpi vec: actual bottom drainage
  Vec :: qflux_subbase_clms   ! seq vec

  Vec :: eflux_subsurf_pfp    ! mpi vec: actual top interface energy flux
  Vec :: eflux_subsurf_clms   ! seq vec
  Vec :: eflux_subbase_pfp    ! mpi vec: actual bottom interface energy flux
  Vec :: eflux_subbase_clms   ! seq vec

  ! net water (with thermal) flow for PFLOTRAN's 3D subsurface domain
  Vec :: qflow_pfp     ! mpi vec (H2O): kgH2O/m3/sec
  Vec :: qflow_clms    ! seq vec
  Vec :: qflowt_pfp    ! mpi vec (H2O)
  Vec :: qflowt_clms   ! seq vec

  ! net heat exchange (no mass) for 3-D subsurface domain
  Vec :: eflow_pfp     ! mpi vec: all non-mass forms of energy exchange rate (MJ/m3/sec)
  Vec :: eflow_clms    ! seq vec



  !---------------------------------------------------------------

  end type clm_pflotran_idata_type

  type(clm_pflotran_idata_type) , public, target , save :: clm_pf_idata
  
  public :: CLMPFLOTRANIDataInit, &
            CLMPFLOTRANIDataCreateVec, &
            CLMPFLOTRANIDataDestroy
  
contains

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataInit()
  ! 
  ! This routine initialized the data transfer type.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none

    nullify(clm_pf_idata%dxclm_global)
    nullify(clm_pf_idata%dyclm_global)
    nullify(clm_pf_idata%dzclm_global)

    clm_pf_idata%npx = 1    !default 'np' for PF mesh decompose is 1x1x1
    clm_pf_idata%npy = 1
    clm_pf_idata%npz = 1
    nullify(clm_pf_idata%clm_lx)
    nullify(clm_pf_idata%clm_ly)
    nullify(clm_pf_idata%clm_lz)

    clm_pf_idata%nzclm_mapped = 0
    clm_pf_idata%nxclm_mapped = 0
    clm_pf_idata%nyclm_mapped = 0

    clm_pf_idata%x0clm_global = 0
    clm_pf_idata%y0clm_global = 0
    clm_pf_idata%z0clm_global = 0

    clm_pf_idata%nlclm_sub = 0
    clm_pf_idata%ngclm_sub = 0
    clm_pf_idata%nlpf_sub  = 0
    clm_pf_idata%ngpf_sub  = 0

    clm_pf_idata%nlclm_2dtop = 0
    clm_pf_idata%ngclm_2dtop = 0
    clm_pf_idata%nlpf_2dtop  = 0
    clm_pf_idata%ngpf_2dtop  = 0

    clm_pf_idata%nlclm_2dbot = 0
    clm_pf_idata%ngclm_2dbot = 0
    clm_pf_idata%nlpf_2dbot  = 0
    clm_pf_idata%ngpf_2dbot  = 0

    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf  = 0
    clm_pf_idata%ngpf_srf  = 0

    !
    clm_pf_idata%zsoil_clmp      = PETSC_NULL_VEC
    clm_pf_idata%zsoil_pfs       = PETSC_NULL_VEC
    clm_pf_idata%zsoil_pfp       = PETSC_NULL_VEC
    clm_pf_idata%zsoil_clms      = PETSC_NULL_VEC
    clm_pf_idata%dxsoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%dxsoil_pfs      = PETSC_NULL_VEC
    clm_pf_idata%dysoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%dysoil_pfs      = PETSC_NULL_VEC
    clm_pf_idata%dzsoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%dzsoil_pfs      = PETSC_NULL_VEC
    clm_pf_idata%xsoil_clmp      = PETSC_NULL_VEC
    clm_pf_idata%xsoil_pfs       = PETSC_NULL_VEC
    clm_pf_idata%ysoil_clmp      = PETSC_NULL_VEC
    clm_pf_idata%ysoil_pfs       = PETSC_NULL_VEC
    clm_pf_idata%zisoil_clmp     = PETSC_NULL_VEC
    clm_pf_idata%zisoil_pfs      = PETSC_NULL_VEC

    clm_pf_idata%area_subsurf_clmp     = PETSC_NULL_VEC
    clm_pf_idata%area_subsurf_pfs      = PETSC_NULL_VEC
    clm_pf_idata%area_subsurf_pfp      = PETSC_NULL_VEC
    clm_pf_idata%area_subsurf_clms     = PETSC_NULL_VEC

    clm_pf_idata%area_top_face_clmp = PETSC_NULL_VEC
    clm_pf_idata%area_top_face_pfs  = PETSC_NULL_VEC
    clm_pf_idata%area_top_face_pfp  = PETSC_NULL_VEC
    clm_pf_idata%area_top_face_clms = PETSC_NULL_VEC

    clm_pf_idata%cellid_clmp     = PETSC_NULL_VEC
    clm_pf_idata%cellid_pfs      = PETSC_NULL_VEC
    clm_pf_idata%cellid_pfp      = PETSC_NULL_VEC
    clm_pf_idata%cellid_clms     = PETSC_NULL_VEC

    clm_pf_idata%cellid_2dtop_clmp     = PETSC_NULL_VEC
    clm_pf_idata%cellid_2dtop_pfs      = PETSC_NULL_VEC
    clm_pf_idata%cellid_2dtop_pfp      = PETSC_NULL_VEC
    clm_pf_idata%cellid_2dtop_clms     = PETSC_NULL_VEC

    !-------------
    clm_pf_idata%head_based = PETSC_TRUE
    clm_pf_idata%pressure_reference = 1.01325d5

    !--------------------------------------------------------------------
    clm_pf_idata%hksat_x_clmp = PETSC_NULL_VEC
    clm_pf_idata%hksat_y_clmp = PETSC_NULL_VEC
    clm_pf_idata%hksat_z_clmp = PETSC_NULL_VEC
    clm_pf_idata%watsat_clmp  = PETSC_NULL_VEC
    clm_pf_idata%watfc_clmp   = PETSC_NULL_VEC
    clm_pf_idata%bulkdensity_dry_clmp = PETSC_NULL_VEC
    clm_pf_idata%effporosity_clmp     = PETSC_NULL_VEC

    clm_pf_idata%tkwet_clmp  = PETSC_NULL_VEC
    clm_pf_idata%tkdry_clmp  = PETSC_NULL_VEC
    clm_pf_idata%tkfrz_clmp  = PETSC_NULL_VEC
    clm_pf_idata%hcvsol_clmp = PETSC_NULL_VEC

    clm_pf_idata%hksat_x_pfs = PETSC_NULL_VEC
    clm_pf_idata%hksat_y_pfs = PETSC_NULL_VEC
    clm_pf_idata%hksat_z_pfs = PETSC_NULL_VEC
    clm_pf_idata%watsat_pfs  = PETSC_NULL_VEC
    clm_pf_idata%watfc_pfs   = PETSC_NULL_VEC
    clm_pf_idata%bulkdensity_dry_pfs = PETSC_NULL_VEC
    clm_pf_idata%effporosity_pfs     = PETSC_NULL_VEC

    clm_pf_idata%tkwet_pfs  = PETSC_NULL_VEC
    clm_pf_idata%tkdry_pfs  = PETSC_NULL_VEC
    clm_pf_idata%tkfrz_pfs  = PETSC_NULL_VEC
    clm_pf_idata%hcvsol_pfs = PETSC_NULL_VEC

    clm_pf_idata%sucsat_clmp = PETSC_NULL_VEC
    clm_pf_idata%bsw_clmp    = PETSC_NULL_VEC
    clm_pf_idata%sucsat_pfs  = PETSC_NULL_VEC
    clm_pf_idata%bsw_pfs     = PETSC_NULL_VEC
   
   !--------------------------------------------------------------------
    clm_pf_idata%sr_pcwmax_pfp   = PETSC_NULL_VEC
    clm_pf_idata%pcwmax_pfp      = PETSC_NULL_VEC
    clm_pf_idata%effporosity_pfp = PETSC_NULL_VEC
    clm_pf_idata%sr_pcwmax_clms  = PETSC_NULL_VEC
    clm_pf_idata%pcwmax_clms     = PETSC_NULL_VEC
    clm_pf_idata%effporosity_clms= PETSC_NULL_VEC

   !--------------------------------------------------------------------
    clm_pf_idata%press_ref_clmp = PETSC_NULL_VEC
    clm_pf_idata%press_ref_pfs  = PETSC_NULL_VEC
    !
    clm_pf_idata%press_clmp    = PETSC_NULL_VEC
    clm_pf_idata%soilpsi_clmp  = PETSC_NULL_VEC
    clm_pf_idata%soillsat_clmp = PETSC_NULL_VEC
    clm_pf_idata%soilisat_clmp = PETSC_NULL_VEC
    clm_pf_idata%soilliq_clmp  = PETSC_NULL_VEC
    clm_pf_idata%soilice_clmp  = PETSC_NULL_VEC
    clm_pf_idata%soilt_clmp    = PETSC_NULL_VEC
    clm_pf_idata%press_pfs      = PETSC_NULL_VEC
    clm_pf_idata%soilpsi_pfs    = PETSC_NULL_VEC
    clm_pf_idata%soillsat_pfs   = PETSC_NULL_VEC
    clm_pf_idata%soilisat_pfs   = PETSC_NULL_VEC
    clm_pf_idata%soilliq_pfs    = PETSC_NULL_VEC
    clm_pf_idata%soilice_pfs    = PETSC_NULL_VEC
    clm_pf_idata%soilt_pfs      = PETSC_NULL_VEC
    !
    clm_pf_idata%press_pfp    = PETSC_NULL_VEC
    clm_pf_idata%soilpsi_pfp  = PETSC_NULL_VEC
    clm_pf_idata%soillsat_pfp = PETSC_NULL_VEC
    clm_pf_idata%soilisat_pfp = PETSC_NULL_VEC
    clm_pf_idata%soilliq_pfp  = PETSC_NULL_VEC
    clm_pf_idata%soilice_pfp  = PETSC_NULL_VEC
    clm_pf_idata%soilt_pfp    = PETSC_NULL_VEC
    clm_pf_idata%press_clms      = PETSC_NULL_VEC
    clm_pf_idata%soilpsi_clms    = PETSC_NULL_VEC
    clm_pf_idata%soillsat_clms   = PETSC_NULL_VEC
    clm_pf_idata%soilisat_clms   = PETSC_NULL_VEC
    clm_pf_idata%soilliq_clms    = PETSC_NULL_VEC
    clm_pf_idata%soilice_clms    = PETSC_NULL_VEC
    clm_pf_idata%soilt_clms      = PETSC_NULL_VEC

    !------------------------------------------------------------------
    clm_pf_idata%N_molecular_weight = 14.0067d0
    clm_pf_idata%C_molecular_weight = 12.0110d0
    nullify(clm_pf_idata%floating_cn_ratio)
    nullify(clm_pf_idata%decomp_element_ratios)
    nullify(clm_pf_idata%ispec_decomp_c)
    nullify(clm_pf_idata%ispec_decomp_n)
    nullify(clm_pf_idata%ispec_decomp_hr)
    nullify(clm_pf_idata%ispec_decomp_nmin)
    nullify(clm_pf_idata%ispec_decomp_nimp)
    nullify(clm_pf_idata%ispec_decomp_nimm)
    nullify(clm_pf_idata%decomp_pool_name)
    nullify(clm_pf_idata%ck_decomp_c)
    nullify(clm_pf_idata%adfactor_ck_c)
    nullify(clm_pf_idata%fr_decomp_c)

    clm_pf_idata%ndecomp_pools    = 0
    clm_pf_idata%ndecomp_elements = 0

    clm_pf_idata%ispec_hrimm   = 0
    clm_pf_idata%ispec_nmin    = 0
    clm_pf_idata%ispec_nimm    = 0
    clm_pf_idata%ispec_nimp    = 0

    clm_pf_idata%ispec_nh4     = 0
    clm_pf_idata%ispec_no3     = 0
    clm_pf_idata%ispec_nh4s    = 0
    clm_pf_idata%ispec_no3s    = 0
    clm_pf_idata%ispec_nh4sorb = 0
    clm_pf_idata%ispec_plantndemand   = 0
    clm_pf_idata%ispec_plantnh4uptake = 0
    clm_pf_idata%ispec_plantno3uptake = 0
    clm_pf_idata%ispec_ngasmin  = 0
    clm_pf_idata%ispec_ngasnitr = 0
    clm_pf_idata%ispec_ngasdeni = 0
    clm_pf_idata%ispec_co2 = 0
    clm_pf_idata%ispec_n2  = 0
    clm_pf_idata%ispec_n2o = 0

   !--------------------------------------------------------------------
    clm_pf_idata%t_scalar_clmp     = PETSC_NULL_VEC
    clm_pf_idata%w_scalar_clmp     = PETSC_NULL_VEC
    clm_pf_idata%o_scalar_clmp     = PETSC_NULL_VEC
    clm_pf_idata%depth_scalar_clmp = PETSC_NULL_VEC
    clm_pf_idata%t_scalar_pfs      = PETSC_NULL_VEC
    clm_pf_idata%w_scalar_pfs      = PETSC_NULL_VEC
    clm_pf_idata%o_scalar_pfs      = PETSC_NULL_VEC
    clm_pf_idata%depth_scalar_pfs  = PETSC_NULL_VEC

    clm_pf_idata%decomp_cpools_vr_clmp = PETSC_NULL_VEC
    clm_pf_idata%decomp_npools_vr_clmp = PETSC_NULL_VEC
    clm_pf_idata%kscalar_decomp_c_clmp = PETSC_NULL_VEC

    clm_pf_idata%smin_no3_vr_clmp      = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4_vr_clmp      = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4sorb_vr_clmp  = PETSC_NULL_VEC

    clm_pf_idata%decomp_cpools_vr_pfs = PETSC_NULL_VEC
    clm_pf_idata%decomp_npools_vr_pfs = PETSC_NULL_VEC
    clm_pf_idata%kscalar_decomp_c_pfs = PETSC_NULL_VEC

    clm_pf_idata%smin_no3_vr_pfs      = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4_vr_pfs      = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4sorb_vr_pfs  = PETSC_NULL_VEC

    clm_pf_idata%rate_decomp_c_clmp         = PETSC_NULL_VEC
    clm_pf_idata%rate_decomp_n_clmp         = PETSC_NULL_VEC
    clm_pf_idata%rate_smin_no3_clmp         = PETSC_NULL_VEC
    clm_pf_idata%rate_smin_nh4_clmp         = PETSC_NULL_VEC
    clm_pf_idata%rate_plantndemand_clmp     = PETSC_NULL_VEC

    clm_pf_idata%rate_decomp_c_pfs         = PETSC_NULL_VEC
    clm_pf_idata%rate_decomp_n_pfs         = PETSC_NULL_VEC
    clm_pf_idata%rate_smin_no3_pfs         = PETSC_NULL_VEC
    clm_pf_idata%rate_smin_nh4_pfs         = PETSC_NULL_VEC
    clm_pf_idata%rate_plantndemand_pfs     = PETSC_NULL_VEC

   !--------------------------------------------------------------------
    ! for C-N states
    clm_pf_idata%decomp_cpools_vr_pfp  = PETSC_NULL_VEC
    clm_pf_idata%decomp_npools_vr_pfp  = PETSC_NULL_VEC
    clm_pf_idata%smin_no3_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4sorb_vr_pfp   = PETSC_NULL_VEC
    clm_pf_idata%decomp_cpools_vr_clms = PETSC_NULL_VEC
    clm_pf_idata%decomp_npools_vr_clms = PETSC_NULL_VEC
    clm_pf_idata%smin_no3_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%smin_nh4sorb_vr_clms  = PETSC_NULL_VEC

    ! for root N extraction calculation
    clm_pf_idata%accextrnh4_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%accextrnh4_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%accextrno3_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%accextrno3_vr_clms      = PETSC_NULL_VEC

    ! for soil hr calculation
    clm_pf_idata%gco2_vr_pfp            = PETSC_NULL_VEC
    clm_pf_idata%gco2_vr_clms           = PETSC_NULL_VEC
    clm_pf_idata%gco2_vr_clmp           = PETSC_NULL_VEC
    clm_pf_idata%gco2_vr_pfs            = PETSC_NULL_VEC

    ! for N2 gas emission calculation
    clm_pf_idata%gn2_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%gn2_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%gn2_vr_clmp      = PETSC_NULL_VEC
    clm_pf_idata%gn2_vr_pfs       = PETSC_NULL_VEC

    ! for N2O gas emission calculation
    clm_pf_idata%gn2o_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%gn2o_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%gn2o_vr_clmp      = PETSC_NULL_VEC
    clm_pf_idata%gn2o_vr_pfs       = PETSC_NULL_VEC

    ! for tracking variables in C-N cycle
    clm_pf_idata%acchr_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%acchr_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%acctothr_vr_pfp    = PETSC_NULL_VEC
    clm_pf_idata%acctothr_vr_clms   = PETSC_NULL_VEC

    clm_pf_idata%accnmin_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%accnmin_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%acctotnmin_vr_pfp    = PETSC_NULL_VEC
    clm_pf_idata%acctotnmin_vr_clms   = PETSC_NULL_VEC

    clm_pf_idata%accnimmp_vr_pfp      = PETSC_NULL_VEC
    clm_pf_idata%accnimmp_vr_clms     = PETSC_NULL_VEC
    clm_pf_idata%acctotnimmp_vr_pfp   = PETSC_NULL_VEC
    clm_pf_idata%acctotnimmp_vr_clms  = PETSC_NULL_VEC

    clm_pf_idata%accnimm_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%accnimm_vr_clms      = PETSC_NULL_VEC
    clm_pf_idata%acctotnimm_vr_pfp    = PETSC_NULL_VEC
    clm_pf_idata%acctotnimm_vr_clms   = PETSC_NULL_VEC

    clm_pf_idata%accngasmin_vr_pfp        = PETSC_NULL_VEC
    clm_pf_idata%accngasmin_vr_clms       = PETSC_NULL_VEC

    clm_pf_idata%accngasnitr_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%accngasnitr_vr_clms      = PETSC_NULL_VEC

    clm_pf_idata%accngasdeni_vr_pfp       = PETSC_NULL_VEC
    clm_pf_idata%accngasdeni_vr_clms      = PETSC_NULL_VEC

    ! aq. chemical species boundary flux

    clm_pf_idata%f_nh4_subsurf_pfp   = PETSC_NULL_VEC
    clm_pf_idata%f_nh4_subsurf_clms  = PETSC_NULL_VEC
    clm_pf_idata%f_nh4_subbase_pfp   = PETSC_NULL_VEC
    clm_pf_idata%f_nh4_subbase_clms  = PETSC_NULL_VEC
    clm_pf_idata%f_no3_subsurf_pfp   = PETSC_NULL_VEC
    clm_pf_idata%f_no3_subsurf_clms  = PETSC_NULL_VEC
    clm_pf_idata%f_no3_subbase_pfp   = PETSC_NULL_VEC
    clm_pf_idata%f_no3_subbase_clms  = PETSC_NULL_VEC

   !--------------------------------------------------------------------

    clm_pf_idata%qflow_clmp = PETSC_NULL_VEC
    clm_pf_idata%qflow_pfs  = PETSC_NULL_VEC
    clm_pf_idata%qflowt_clmp= PETSC_NULL_VEC
    clm_pf_idata%qflowt_pfs = PETSC_NULL_VEC
    clm_pf_idata%eflow_clmp = PETSC_NULL_VEC
    clm_pf_idata%eflow_pfs  = PETSC_NULL_VEC

    clm_pf_idata%press_maxponding_clmp = PETSC_NULL_VEC
    clm_pf_idata%press_maxponding_pfs  = PETSC_NULL_VEC
    clm_pf_idata%press_subsurf_clmp = PETSC_NULL_VEC
    clm_pf_idata%press_subsurf_pfs  = PETSC_NULL_VEC
    clm_pf_idata%press_subbase_clmp = PETSC_NULL_VEC
    clm_pf_idata%press_subbase_pfs  = PETSC_NULL_VEC

    clm_pf_idata%qfluxw_subsurf_clmp  = PETSC_NULL_VEC
    clm_pf_idata%qfluxw_subsurf_pfs   = PETSC_NULL_VEC
    clm_pf_idata%qfluxev_subsurf_clmp = PETSC_NULL_VEC
    clm_pf_idata%qfluxev_subsurf_pfs  = PETSC_NULL_VEC
    clm_pf_idata%qfluxw_subbase_clmp  = PETSC_NULL_VEC
    clm_pf_idata%qfluxw_subbase_pfs   = PETSC_NULL_VEC

    clm_pf_idata%gtemp_subsurf_clmp = PETSC_NULL_VEC
    clm_pf_idata%gtemp_subsurf_pfs  = PETSC_NULL_VEC
    clm_pf_idata%eflux_subsurf_clmp = PETSC_NULL_VEC
    clm_pf_idata%eflux_subsurf_pfs  = PETSC_NULL_VEC
    clm_pf_idata%efluxr_subsurf_clmp= PETSC_NULL_VEC
    clm_pf_idata%efluxr_subsurf_pfs = PETSC_NULL_VEC
    clm_pf_idata%efluxl_subsurf_clmp= PETSC_NULL_VEC
    clm_pf_idata%efluxl_subsurf_pfs = PETSC_NULL_VEC
    clm_pf_idata%gtemp_subbase_clmp = PETSC_NULL_VEC
    clm_pf_idata%gtemp_subbase_pfs  = PETSC_NULL_VEC
    clm_pf_idata%eflux_subbase_clmp = PETSC_NULL_VEC
    clm_pf_idata%eflux_subbase_pfs  = PETSC_NULL_VEC

    !---------------
    ! water/energy boundary flux
    clm_pf_idata%qinfl_subsurf_pfp   = PETSC_NULL_VEC
    clm_pf_idata%qinfl_subsurf_clms  = PETSC_NULL_VEC
    clm_pf_idata%qsurf_subsurf_pfp   = PETSC_NULL_VEC
    clm_pf_idata%qsurf_subsurf_clms  = PETSC_NULL_VEC
    clm_pf_idata%qflux_subbase_pfp   = PETSC_NULL_VEC
    clm_pf_idata%qflux_subbase_clms  = PETSC_NULL_VEC
    clm_pf_idata%eflux_subsurf_pfp   = PETSC_NULL_VEC
    clm_pf_idata%eflux_subsurf_clms  = PETSC_NULL_VEC
    clm_pf_idata%eflux_subbase_pfp   = PETSC_NULL_VEC
    clm_pf_idata%eflux_subbase_clms  = PETSC_NULL_VEC

    clm_pf_idata%qflow_pfp   = PETSC_NULL_VEC
    clm_pf_idata%qflow_clms  = PETSC_NULL_VEC
    clm_pf_idata%qflowt_pfp  = PETSC_NULL_VEC
    clm_pf_idata%qflowt_clms = PETSC_NULL_VEC
    clm_pf_idata%eflow_pfp   = PETSC_NULL_VEC
    clm_pf_idata%eflow_clms  = PETSC_NULL_VEC

  end subroutine CLMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataCreateVec(mycomm)
  ! 
  ! This routine creates PETSc vectors required for data transfer between
  ! CLM and PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank

    call MPI_Comm_rank(mycomm,rank, ierr)

    ! The following block of data definition is for THC coupled clm-pflotran (Currently ONLY subsurface or soil)
    !
    !NOTES (fmy): From mpi vecs To seq. vecs for passing data IS in one-way only at this momment
    !             (1) First, here will create 4 sets of 3D/2D vecs.
    !             (2) then, below will copy these vecs to create vecs for other variables.

    ! -------- FOR CLM (mpi) ==> PFLOTRAN (seq)
    ! CLM(mpi)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%zsoil_clmp,ierr)             ! 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%zsoil_clmp,0.d0,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_2dtop,PETSC_DECIDE,clm_pf_idata%area_subsurf_clmp,ierr)     ! 2D top-cells of 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%area_subsurf_clmp,0.d0,ierr)
    ! PFLOTRAN(seq)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%zsoil_pfs,ierr)                   ! 3D Subsurface CLM
    call VecSet(clm_pf_idata%zsoil_pfs,0.d0,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_2dtop,clm_pf_idata%area_subsurf_pfs,ierr)           ! 2D top-cells of 3D Subsurface CLM
    call VecSet(clm_pf_idata%area_subsurf_pfs,0.d0,ierr)

    ! -------- FOR PFLOTRAN (mpi) ==> CLM (seq)
    ! PFLOTRAN(mpi)
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%zsoil_pfp,ierr)               ! 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%zsoil_pfp,0.d0,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_2dtop,PETSC_DECIDE,clm_pf_idata%area_subsurf_pfp,ierr)     ! 2D top-cells of 3D Subsurface PFLOTRAN
    call VecSet(clm_pf_idata%area_subsurf_pfp,0.d0,ierr)
    ! CLM(seq)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%zsoil_clms,ierr)                 ! 3D Subsurface CLM
    call VecSet(clm_pf_idata%zsoil_clms,0.d0,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_2dtop,clm_pf_idata%area_subsurf_clms,ierr)       ! 2D top-cells of 3D Subsurface CLM
    call VecSet(clm_pf_idata%area_subsurf_clms,0.d0,ierr)



    !
    ! ----------------------------------------------------------------------------
    !
    ! I. CONSTANTS, INITIALS
    !
    !-----------------------------CLM ==> PFLOTRAN
    ! (by copying) Create MPI Vectors for CLM ---------------------------------
    ! soil dimensions
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%xsoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%ysoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%zisoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%dxsoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%dysoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%dzsoil_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%area_top_face_clmp,ierr)

    ! soil cell ids (3D) / surface cell ids (2D)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%cellid_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%cellid_2dtop_clmp,ierr)

    ! soil physical properties (3D)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hksat_x_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hksat_y_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hksat_z_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%watsat_clmp,ierr)       ! total vwc at saturation (total 'porosity')
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%watfc_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%bulkdensity_dry_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%effporosity_clmp,ierr)     ! this may/may not same as 'bd'/'watsat' above

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%tkwet_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%tkdry_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%tkfrz_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%hcvsol_clmp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%sucsat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%bsw_clmp,ierr)

    ! TH states (3D)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%press_ref_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%press_clmp,ierr)        ! this depends upon 'reference pressure'
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%soilpsi_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%soillsat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%soilisat_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%soilliq_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%soilice_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%soilt_clmp,ierr)

    ! (by copying) Create Seq. Vectors for PFLOTRAN  ----------------------
    ! soil dimensions
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%xsoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%ysoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%zisoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%dxsoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%dysoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%dzsoil_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%area_top_face_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%cellid_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%cellid_2dtop_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hksat_x_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hksat_y_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hksat_z_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%watsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%watfc_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%bulkdensity_dry_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%effporosity_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%tkwet_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%tkdry_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%tkfrz_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%hcvsol_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%sucsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%bsw_pfs,ierr)

     ! TH states (3D)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%press_ref_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%press_pfs,ierr)        ! this depends upon 'reference pressure'
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%soilpsi_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%soillsat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%soilisat_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%soilliq_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%soilice_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%soilt_pfs,ierr)

    !-----------------------------PFLOTRAN ==> CLM
    ! (by copying) Create MPI Vectors for PFLOTRAN ------------
    ! soil dimensions
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%area_top_face_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%cellid_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%cellid_2dtop_pfp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%sr_pcwmax_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%pcwmax_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%effporosity_pfp,ierr)

     ! TH states (3D)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%press_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%soilpsi_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%soillsat_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%soilisat_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%soilliq_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%soilice_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%soilt_pfp,ierr)

    ! (by copying) create Seq. Vectors for CLM ---------
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%area_top_face_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%cellid_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%cellid_2dtop_clms,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%sr_pcwmax_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%pcwmax_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%effporosity_clms,ierr)

     ! TH states (3D)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%press_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%soilpsi_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%soillsat_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%soilisat_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%soilliq_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%soilice_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%soilt_clms,ierr)

    !--------------------------------------------------------------------------------------------------------------------------
    !
    ! II. BGC VARIABLES
    !
    ! BGC state variables: 3D subsurface CLM ---to--- 3D subsurface PFLOTRAN (e.g., initialization or restarting)
    !-----------------------------CLM ==> PFLOTRAN
    ! MPI Vecs for CLM
    call VecCreateMPI(mycomm, clm_pf_idata%ndecomp_pools*clm_pf_idata%nlclm_sub,   &    ! no. of decomp_pools X 3D Subsurface cells
         PETSC_DECIDE,clm_pf_idata%decomp_cpools_vr_clmp,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_clmp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clmp, clm_pf_idata%decomp_npools_vr_clmp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%kscalar_decomp_c_clmp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%t_scalar_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%w_scalar_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%o_scalar_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%depth_scalar_clmp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%smin_no3_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%smin_nh4_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%smin_nh4sorb_vr_clmp,ierr)

    ! Seq. Vecs for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ndecomp_pools*clm_pf_idata%ngpf_sub, &   ! no. of decomp_pools X 3D Subsurface cells
          clm_pf_idata%decomp_cpools_vr_pfs,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_pfs,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfs, clm_pf_idata%decomp_npools_vr_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs, clm_pf_idata%kscalar_decomp_c_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%t_scalar_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%w_scalar_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%o_scalar_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%depth_scalar_pfs,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%smin_no3_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%smin_nh4_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%smin_nh4sorb_vr_pfs,ierr)

    ! BGC/TH interface source/sink (rate): 3D subsurface CLM ---to--- 3D subsurface PFLOTRAN
    ! MPI Vecs for CLM
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clmp,clm_pf_idata%rate_decomp_c_clmp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clmp,clm_pf_idata%rate_decomp_n_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%rate_smin_no3_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%rate_smin_nh4_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%rate_plantndemand_clmp,ierr)
    ! Seq. Vecs for PFLOTRAN
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfs,clm_pf_idata%rate_decomp_c_pfs,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfs,clm_pf_idata%rate_decomp_n_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%rate_smin_no3_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%rate_smin_nh4_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%rate_plantndemand_pfs,ierr)

    ! MPI Vecs for CLM to pass reset aq. conc back to PF
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%gco2_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%gn2_vr_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%gn2o_vr_clmp,ierr)
    ! Seq. Vecs for PFLOTRAN to get reset aq. conc back from CLM
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%gco2_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%gn2_vr_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%gn2o_vr_pfs,ierr)

    !-----------------------------PFLOTRAN ==> CLM
    ! BGC state variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecCreateMPI(mycomm, clm_pf_idata%ndecomp_pools*clm_pf_idata%nlpf_sub,   &    ! no. of decomp_pools X 3D Subsurface cells
         PETSC_DECIDE,clm_pf_idata%decomp_cpools_vr_pfp,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_pfp,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfp, clm_pf_idata%decomp_npools_vr_pfp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%smin_no3_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%smin_nh4_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%smin_nh4sorb_vr_pfp,ierr)

    ! Seq. Vecs for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ndecomp_pools*clm_pf_idata%ngclm_sub, &   ! no. of decomp_pools X 3D Subsurface cells
          clm_pf_idata%decomp_cpools_vr_clms,ierr)
    call VecSet(clm_pf_idata%decomp_cpools_vr_clms,0.d0,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clms,clm_pf_idata%decomp_npools_vr_clms,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%smin_no3_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%smin_nh4_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%smin_nh4sorb_vr_clms,ierr)


    ! BGC flux variables: 3D subsurface PFLOTRAN ---to--- 3D subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%gco2_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%gn2_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%gn2o_vr_pfp,ierr)
    !
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%accextrnh4_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%accextrno3_vr_pfp,ierr)
    !
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfp,clm_pf_idata%acchr_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfp,clm_pf_idata%accnmin_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfp,clm_pf_idata%accnimmp_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_pfp,clm_pf_idata%accnimm_vr_pfp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%acctothr_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%acctotnmin_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%acctotnimmp_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%acctotnimm_vr_pfp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%accngasmin_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%accngasnitr_vr_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%accngasdeni_vr_pfp,ierr)

    ! Seq. Vecs for CLM
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%gco2_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%gn2_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%gn2o_vr_clms,ierr)
    !
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%accextrnh4_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%accextrno3_vr_clms,ierr)
    !
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clms,clm_pf_idata%acchr_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clms,clm_pf_idata%accnmin_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clms,clm_pf_idata%accnimmp_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%decomp_cpools_vr_clms,clm_pf_idata%accnimm_vr_clms,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%acctothr_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%acctotnmin_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%acctotnimmp_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%acctotnimm_vr_clms,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%accngasmin_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%accngasnitr_vr_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%accngasdeni_vr_clms,ierr)

    ! BC flow variables: 2D faces of subsurface PFLOTRAN ---to--- 2D faces of subsurface CLM
    ! MPI Vecs for PFLOTRAN
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_nh4_subsurf_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_no3_subsurf_pfp,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_nh4_subbase_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%f_no3_subbase_pfp,ierr)
    ! Seq. Vecs for CLM
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_nh4_subsurf_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_no3_subsurf_clms,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_nh4_subbase_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%f_no3_subbase_clms,ierr)

    !
    !--------------------------------------------------------------------------------------------------------------------------
    ! THERMAL-HYDROLOGY VARIABLES
    !
    ! --------- For TH data transfer from CLM to PFLOTRAN
    ! (by copying) Create mpi Vectors for CLM  ----------------------
    ! TH Src/Sink (3D)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%qflow_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%qflowt_clmp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clmp,clm_pf_idata%eflow_clmp,ierr)

    ! TH top BC (2D)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%gtemp_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%qfluxw_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%qfluxev_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%eflux_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%efluxr_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%efluxl_subsurf_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_maxponding_clmp,ierr)

    ! TH bottom BC (2D): supposing SAME bottom-cell numbers as top-cells!
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%press_subbase_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%gtemp_subbase_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%qfluxw_subbase_clmp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clmp,clm_pf_idata%eflux_subbase_clmp,ierr)

    ! (by copying) Create Seq. Vectors for PFLOTRAN  ----------------------
    ! TH Src/Sink (3D)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%qflow_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%qflowt_pfs,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfs,clm_pf_idata%eflow_pfs,ierr)

    ! TH top BC (2D)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_subsurf_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%gtemp_subsurf_pfs,ierr)            ! T
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%qfluxw_subsurf_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%qfluxev_subsurf_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%eflux_subsurf_pfs,ierr)            ! T
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%efluxr_subsurf_pfs,ierr)            ! T
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%efluxl_subsurf_pfs,ierr)            ! T
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_maxponding_pfs,ierr)

    ! TH bottom BC (2D): supposing SAME bot-cell numbers as top-cells!
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%press_subbase_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%gtemp_subbase_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%qfluxw_subbase_pfs,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfs,clm_pf_idata%eflux_subbase_pfs,ierr)            ! T

    ! --------- For TH data transfer from PFLOTRAN to CLM
    ! (by copying) Create MPI Vectors for PFLOTRAN  ----------------------
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%qinfl_subsurf_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%qsurf_subsurf_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%eflux_subsurf_pfp,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%qflux_subbase_pfp,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_pfp,clm_pf_idata%eflux_subbase_pfp,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%qflow_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%qflowt_pfp,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_pfp,clm_pf_idata%eflow_pfp,ierr)

    ! (by copying) Create Seq. Vectors for CLM  ----------------------
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%qinfl_subsurf_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%qsurf_subsurf_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%eflux_subsurf_clms,ierr)
    ! the following assumes SAME bot-cell number as top-cells
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%qflux_subbase_clms,ierr)
    call VecDuplicate(clm_pf_idata%area_subsurf_clms,clm_pf_idata%eflux_subbase_clms,ierr)

    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%qflow_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%qflowt_clms,ierr)
    call VecDuplicate(clm_pf_idata%zsoil_clms,clm_pf_idata%eflow_clms,ierr)

    !--------------------------------------------------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataDestroy()
  ! 
  ! This routine destroys PETSc vectors that were created for data transfer.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! Revised by Fengming Yuan, CCSI-ORNL @May-2015
  
    implicit none
    
    PetscErrorCode :: ierr

    if (associated(clm_pf_idata%dxclm_global)) &
    deallocate(clm_pf_idata%dxclm_global)
    if (associated(clm_pf_idata%dyclm_global)) &
    deallocate(clm_pf_idata%dyclm_global)
    if (associated(clm_pf_idata%dzclm_global)) &
    deallocate(clm_pf_idata%dzclm_global)

    !----------------------------------------------------------------------------------

    if(clm_pf_idata%zsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_clmp,ierr)
    if(clm_pf_idata%zsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_pfs,ierr)
    if(clm_pf_idata%zsoil_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_pfp,ierr)
    if(clm_pf_idata%zsoil_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zsoil_clms,ierr)

    if(clm_pf_idata%xsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%xsoil_clmp,ierr)
    if(clm_pf_idata%xsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%xsoil_pfs,ierr)
    if(clm_pf_idata%ysoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%ysoil_clmp,ierr)
    if(clm_pf_idata%ysoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%ysoil_pfs,ierr)
    if(clm_pf_idata%zisoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zisoil_clmp,ierr)
    if(clm_pf_idata%zisoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%zisoil_pfs,ierr)

    if(clm_pf_idata%dxsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dxsoil_clmp,ierr)
    if(clm_pf_idata%dxsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dxsoil_pfs,ierr)
    if(clm_pf_idata%dysoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dysoil_clmp,ierr)
    if(clm_pf_idata%dysoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dysoil_pfs,ierr)
    if(clm_pf_idata%dzsoil_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dzsoil_clmp,ierr)
    if(clm_pf_idata%dzsoil_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%dzsoil_pfs,ierr)

    if(clm_pf_idata%area_subsurf_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_clmp,ierr)
    if(clm_pf_idata%area_subsurf_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_pfs,ierr)
    if(clm_pf_idata%area_subsurf_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_pfp,ierr)
    if(clm_pf_idata%area_subsurf_clms  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_subsurf_clms,ierr)

    if(clm_pf_idata%area_top_face_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_clmp,ierr)
    if(clm_pf_idata%area_top_face_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_pfs,ierr)
    if(clm_pf_idata%area_top_face_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_pfp,ierr)
    if(clm_pf_idata%area_top_face_clms  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%area_top_face_clms,ierr)

    !----
    if(clm_pf_idata%cellid_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_clmp,ierr)
    if(clm_pf_idata%cellid_2dtop_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_clmp,ierr)
    if(clm_pf_idata%cellid_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_pfs,ierr)
    if(clm_pf_idata%cellid_2dtop_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_pfs,ierr)
    if(clm_pf_idata%hksat_x_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_x_clmp,ierr)
    if(clm_pf_idata%hksat_y_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_y_clmp,ierr)
    if(clm_pf_idata%hksat_z_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_z_clmp,ierr)
    if(clm_pf_idata%sucsat_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%sucsat_clmp,ierr)
    if(clm_pf_idata%watsat_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watsat_clmp,ierr)
    if(clm_pf_idata%bsw_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bsw_clmp,ierr)
    if(clm_pf_idata%watfc_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watfc_clmp,ierr)
    if(clm_pf_idata%bulkdensity_dry_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bulkdensity_dry_clmp,ierr)

    if(clm_pf_idata%cellid_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_pfp,ierr)
    if(clm_pf_idata%cellid_2dtop_pfp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_pfp,ierr)
    if(clm_pf_idata%cellid_clms  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_clms,ierr)
    if(clm_pf_idata%cellid_2dtop_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%cellid_2dtop_clms,ierr)
    if(clm_pf_idata%tkwet_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkwet_clmp,ierr)
    if(clm_pf_idata%tkdry_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkdry_clmp,ierr)
    if(clm_pf_idata%tkfrz_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkfrz_clmp,ierr)
    if(clm_pf_idata%hcvsol_clmp  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hcvsol_clmp,ierr)

    if(clm_pf_idata%hksat_x_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_x_pfs,ierr)
    if(clm_pf_idata%hksat_y_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_y_pfs,ierr)
    if(clm_pf_idata%hksat_z_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hksat_z_pfs,ierr)
    if(clm_pf_idata%sucsat_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%sucsat_pfs,ierr)
    if(clm_pf_idata%watsat_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watsat_pfs,ierr)
    if(clm_pf_idata%bsw_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bsw_pfs,ierr)
    if(clm_pf_idata%watfc_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%watfc_pfs,ierr)
    if(clm_pf_idata%bulkdensity_dry_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%bulkdensity_dry_pfs,ierr)
    if(clm_pf_idata%effporosity_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_clmp,ierr)
    if(clm_pf_idata%effporosity_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_pfs,ierr)

    !----
    if(clm_pf_idata%tkwet_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkwet_pfs,ierr)
    if(clm_pf_idata%tkdry_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkdry_pfs,ierr)
    if(clm_pf_idata%tkfrz_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%tkfrz_pfs,ierr)
    if(clm_pf_idata%hcvsol_pfs  /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%hcvsol_pfs,ierr)

    ! -----
    if(clm_pf_idata%sr_pcwmax_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%sr_pcwmax_pfp,ierr)
    if(clm_pf_idata%pcwmax_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%pcwmax_pfp,ierr)
    if(clm_pf_idata%sr_pcwmax_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%sr_pcwmax_clms,ierr)
    if(clm_pf_idata%pcwmax_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%pcwmax_clms,ierr)
    if(clm_pf_idata%effporosity_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_pfp,ierr)
    if(clm_pf_idata%effporosity_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%effporosity_clms,ierr)

    !----
    if(clm_pf_idata%press_ref_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_ref_clmp,ierr)
    if(clm_pf_idata%press_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_clmp,ierr)
    if(clm_pf_idata%soilpsi_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilpsi_clmp,ierr)
    if(clm_pf_idata%soillsat_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soillsat_clmp,ierr)
    if(clm_pf_idata%soilisat_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilisat_clmp,ierr)
    if(clm_pf_idata%soilliq_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilliq_clmp,ierr)
    if(clm_pf_idata%soilice_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilice_clmp,ierr)
    if(clm_pf_idata%soilt_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilt_clmp,ierr)
    
    if(clm_pf_idata%press_ref_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%press_ref_pfs,ierr)
    if(clm_pf_idata%press_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%press_pfs,ierr)
    if(clm_pf_idata%soilpsi_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilpsi_pfs,ierr)
    if(clm_pf_idata%soillsat_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soillsat_pfs,ierr)
    if(clm_pf_idata%soilisat_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilisat_pfs,ierr)
    if(clm_pf_idata%soilliq_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilliq_pfs,ierr)
    if(clm_pf_idata%soilice_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilice_pfs,ierr)
    if(clm_pf_idata%soilt_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilt_pfs,ierr)

    !----
    if(clm_pf_idata%press_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_pfp,ierr)
    if(clm_pf_idata%soilpsi_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilpsi_pfp,ierr)
    if(clm_pf_idata%soillsat_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soillsat_pfp,ierr)
    if(clm_pf_idata%soilisat_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilisat_pfp,ierr)
    if(clm_pf_idata%soilliq_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilliq_pfp,ierr)
    if(clm_pf_idata%soilice_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilice_pfp,ierr)
    if(clm_pf_idata%soilt_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%soilt_pfp,ierr)

    if(clm_pf_idata%press_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%press_clms,ierr)
    if(clm_pf_idata%soilpsi_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilpsi_clms,ierr)
    if(clm_pf_idata%soillsat_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soillsat_clms,ierr)
    if(clm_pf_idata%soilisat_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilisat_clms,ierr)
    if(clm_pf_idata%soilliq_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilliq_clms,ierr)
    if(clm_pf_idata%soilice_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilice_clms,ierr)
    if(clm_pf_idata%soilt_clms /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%soilt_clms,ierr)

    ! -----------------------------------------------------------------------------------------------------------

    if (associated(clm_pf_idata%decomp_element_ratios)) &
    deallocate(clm_pf_idata%decomp_element_ratios)
    if (associated(clm_pf_idata%floating_cn_ratio)) &
    deallocate(clm_pf_idata%floating_cn_ratio)
    if (associated(clm_pf_idata%ck_decomp_c)) &
    deallocate(clm_pf_idata%ck_decomp_c)
    if (associated(clm_pf_idata%adfactor_ck_c)) &
    deallocate(clm_pf_idata%adfactor_ck_c)
    if (associated(clm_pf_idata%fr_decomp_c)) &
    deallocate(clm_pf_idata%fr_decomp_c)

    if (associated(clm_pf_idata%ispec_decomp_c)) &
    deallocate(clm_pf_idata%ispec_decomp_c)
    if (associated(clm_pf_idata%ispec_decomp_n)) &
    deallocate(clm_pf_idata%ispec_decomp_n)
    if (associated(clm_pf_idata%ispec_decomp_hr)) &
    deallocate(clm_pf_idata%ispec_decomp_hr)
    if (associated(clm_pf_idata%ispec_decomp_nmin)) &
    deallocate(clm_pf_idata%ispec_decomp_nmin)
    if (associated(clm_pf_idata%ispec_decomp_nimp)) &
    deallocate(clm_pf_idata%ispec_decomp_nimp)
    if (associated(clm_pf_idata%ispec_decomp_nimm)) &
    deallocate(clm_pf_idata%ispec_decomp_nimm)
    if (associated(clm_pf_idata%decomp_pool_name)) &
    deallocate(clm_pf_idata%decomp_pool_name)

    ! soil C/N pools (initial)
    if(clm_pf_idata%decomp_cpools_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_clmp,ierr)
    if(clm_pf_idata%decomp_npools_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_clmp,ierr)

    if(clm_pf_idata%kscalar_decomp_c_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%kscalar_decomp_c_clmp,ierr)

    if(clm_pf_idata%t_scalar_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%t_scalar_clmp,ierr)
    if(clm_pf_idata%w_scalar_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%w_scalar_clmp,ierr)
    if(clm_pf_idata%o_scalar_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%o_scalar_clmp,ierr)
    if(clm_pf_idata%depth_scalar_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%depth_scalar_clmp,ierr)

    if(clm_pf_idata%smin_no3_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clmp,ierr)
    if(clm_pf_idata%smin_nh4_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clmp,ierr)
    if(clm_pf_idata%smin_nh4sorb_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_nh4sorb_vr_clmp,ierr)

    !
    if(clm_pf_idata%decomp_cpools_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_pfs,ierr)
    if(clm_pf_idata%decomp_npools_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_pfs,ierr)

    if(clm_pf_idata%kscalar_decomp_c_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%kscalar_decomp_c_pfs,ierr)

    if(clm_pf_idata%t_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%t_scalar_pfs,ierr)
    if(clm_pf_idata%w_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%w_scalar_pfs,ierr)
    if(clm_pf_idata%o_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%o_scalar_pfs,ierr)
    if(clm_pf_idata%depth_scalar_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%depth_scalar_pfs,ierr)

    if(clm_pf_idata%smin_no3_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfs,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfs,ierr)
    if(clm_pf_idata%smin_nh4sorb_vr_pfs /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%smin_nh4sorb_vr_pfs,ierr)


    ! soil C/N fluxes at interface (source/sink)
    if(clm_pf_idata%rate_decomp_c_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_decomp_c_clmp,ierr)
    if(clm_pf_idata%rate_decomp_n_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_decomp_n_clmp,ierr)

    if(clm_pf_idata%rate_plantndemand_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_plantndemand_clmp,ierr)
    if(clm_pf_idata%rate_smin_no3_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_smin_no3_clmp,ierr)
    if(clm_pf_idata%rate_smin_nh4_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_smin_nh4_clmp,ierr)

    if(clm_pf_idata%rate_decomp_c_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_decomp_c_pfs,ierr)
    if(clm_pf_idata%rate_decomp_n_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_decomp_n_pfs,ierr)

    if(clm_pf_idata%rate_plantndemand_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_plantndemand_pfs,ierr)
    if(clm_pf_idata%rate_smin_no3_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_smin_no3_pfs,ierr)
    if(clm_pf_idata%rate_smin_nh4_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%rate_smin_nh4_pfs,ierr)

    !------
    if(clm_pf_idata%decomp_cpools_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_pfp,ierr)
    if(clm_pf_idata%decomp_npools_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_pfp,ierr)
    if(clm_pf_idata%smin_no3_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_pfp,ierr)
    if(clm_pf_idata%smin_nh4_vr_pfp /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%smin_nh4_vr_pfp,ierr)
    if(clm_pf_idata%smin_nh4sorb_vr_pfp /= PETSC_NULL_VEC) &
      call VecDestroy(clm_pf_idata%smin_nh4sorb_vr_pfp,ierr)

    if(clm_pf_idata%decomp_cpools_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_cpools_vr_clms,ierr)
    if(clm_pf_idata%decomp_npools_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%decomp_npools_vr_clms,ierr)
    if(clm_pf_idata%smin_no3_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_no3_vr_clms,ierr)
    if(clm_pf_idata%smin_nh4_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_nh4_vr_clms,ierr)
    if(clm_pf_idata%smin_nh4sorb_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%smin_nh4sorb_vr_clms,ierr)

    ! -----------
    if(clm_pf_idata%accextrnh4_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accextrnh4_vr_pfp,ierr)
    if(clm_pf_idata%accextrno3_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accextrno3_vr_pfp,ierr)
    if(clm_pf_idata%accextrnh4_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accextrnh4_vr_clms,ierr)
    if(clm_pf_idata%accextrno3_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accextrno3_vr_clms,ierr)

    if(clm_pf_idata%gco2_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gco2_vr_pfp,ierr)
    if(clm_pf_idata%gco2_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gco2_vr_clms,ierr)
    if(clm_pf_idata%gco2_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gco2_vr_clmp,ierr)
    if(clm_pf_idata%gco2_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gco2_vr_pfs,ierr)

    if(clm_pf_idata%gn2_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2_vr_pfp,ierr)
    if(clm_pf_idata%gn2_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2_vr_clms,ierr)
    if(clm_pf_idata%gn2_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2_vr_clmp,ierr)
    if(clm_pf_idata%gn2_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2_vr_pfs,ierr)

    if(clm_pf_idata%gn2o_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2o_vr_pfp,ierr)
    if(clm_pf_idata%gn2o_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2o_vr_clms,ierr)
    if(clm_pf_idata%gn2o_vr_clmp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2o_vr_clmp,ierr)
    if(clm_pf_idata%gn2o_vr_pfs /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gn2o_vr_pfs,ierr)

    if(clm_pf_idata%acchr_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acchr_vr_pfp,ierr)
    if(clm_pf_idata%acchr_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acchr_vr_clms,ierr)

    if(clm_pf_idata%acctothr_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctothr_vr_pfp,ierr)
    if(clm_pf_idata%acctothr_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctothr_vr_clms,ierr)

    if(clm_pf_idata%accnmin_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accnmin_vr_pfp,ierr)
    if(clm_pf_idata%accnmin_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accnmin_vr_clms,ierr)

    if(clm_pf_idata%acctotnmin_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctotnmin_vr_pfp,ierr)
    if(clm_pf_idata%acctotnmin_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctotnmin_vr_clms,ierr)

    if(clm_pf_idata%accnimmp_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accnimmp_vr_pfp,ierr)
    if(clm_pf_idata%accnimmp_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accnimmp_vr_clms,ierr)

    if(clm_pf_idata%acctotnimmp_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctotnimmp_vr_pfp,ierr)
    if(clm_pf_idata%acctotnimmp_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctotnimmp_vr_clms,ierr)

    if(clm_pf_idata%accnimm_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accnimm_vr_pfp,ierr)
    if(clm_pf_idata%accnimm_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accnimm_vr_clms,ierr)

    if(clm_pf_idata%acctotnimm_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctotnimm_vr_pfp,ierr)
    if(clm_pf_idata%acctotnimm_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%acctotnimm_vr_clms,ierr)

    if(clm_pf_idata%accngasmin_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accngasmin_vr_pfp,ierr)
    if(clm_pf_idata%accngasmin_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accngasmin_vr_clms,ierr)

    if(clm_pf_idata%accngasnitr_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accngasnitr_vr_pfp,ierr)
    if(clm_pf_idata%accngasnitr_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accngasnitr_vr_clms,ierr)

    if(clm_pf_idata%accngasdeni_vr_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accngasdeni_vr_pfp,ierr)
    if(clm_pf_idata%accngasdeni_vr_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%accngasdeni_vr_clms,ierr)

    !-------
    if(clm_pf_idata%f_nh4_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_nh4_subsurf_pfp,ierr)
    if(clm_pf_idata%f_nh4_subsurf_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_nh4_subsurf_clms,ierr)
    if(clm_pf_idata%f_nh4_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_nh4_subbase_pfp,ierr)
    if(clm_pf_idata%f_nh4_subbase_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_nh4_subbase_clms,ierr)

    if(clm_pf_idata%f_no3_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_no3_subsurf_pfp,ierr)
    if(clm_pf_idata%f_no3_subsurf_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_no3_subsurf_clms,ierr)
    if(clm_pf_idata%f_no3_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_no3_subbase_pfp,ierr)
    if(clm_pf_idata%f_no3_subbase_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%f_no3_subbase_clms,ierr)
    !
    ! -----------------------------------------------------------------------------------------------------------
    !-----
    if(clm_pf_idata%qflow_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflow_clmp,ierr)
    if(clm_pf_idata%qflow_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflow_pfs,ierr)
    if(clm_pf_idata%qflowt_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflowt_clmp,ierr)
    if(clm_pf_idata%qflowt_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflowt_pfs,ierr)
    if(clm_pf_idata%eflow_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflow_clmp,ierr)
    if(clm_pf_idata%eflow_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflow_pfs,ierr)

    !-----
    if(clm_pf_idata%press_maxponding_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_maxponding_clmp,ierr)
    if(clm_pf_idata%press_maxponding_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_maxponding_pfs,ierr)
    if(clm_pf_idata%press_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_subsurf_clmp,ierr)
    if(clm_pf_idata%press_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_subsurf_pfs,ierr)
    if(clm_pf_idata%press_subbase_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_subbase_clmp,ierr)
    if(clm_pf_idata%press_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%press_subbase_pfs,ierr)
    if(clm_pf_idata%qfluxw_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qfluxw_subsurf_clmp,ierr)
    if(clm_pf_idata%qfluxw_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qfluxw_subsurf_pfs,ierr)
    if(clm_pf_idata%qfluxw_subbase_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qfluxw_subbase_clmp,ierr)
    if(clm_pf_idata%qfluxw_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qfluxw_subbase_pfs,ierr)
    if(clm_pf_idata%qfluxev_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qfluxev_subsurf_clmp,ierr)
    if(clm_pf_idata%qfluxev_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qfluxev_subsurf_pfs,ierr)

    if(clm_pf_idata%efluxr_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%efluxr_subsurf_clmp,ierr)
    if(clm_pf_idata%efluxr_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%efluxr_subsurf_pfs,ierr)
    if(clm_pf_idata%efluxl_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%efluxl_subsurf_clmp,ierr)
    if(clm_pf_idata%efluxl_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%efluxl_subsurf_pfs,ierr)
    if(clm_pf_idata%eflux_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subsurf_clmp,ierr)
    if(clm_pf_idata%eflux_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subsurf_pfs,ierr)
    if(clm_pf_idata%gtemp_subsurf_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gtemp_subsurf_clmp,ierr)
    if(clm_pf_idata%gtemp_subsurf_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gtemp_subsurf_pfs,ierr)
    if(clm_pf_idata%eflux_subbase_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subbase_clmp,ierr)
    if(clm_pf_idata%eflux_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subbase_pfs,ierr)
    if(clm_pf_idata%gtemp_subbase_clmp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gtemp_subbase_clmp,ierr)
    if(clm_pf_idata%gtemp_subbase_pfs  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%gtemp_subbase_pfs,ierr)

    !------------------
    if(clm_pf_idata%qinfl_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qinfl_subsurf_pfp,ierr)
    if(clm_pf_idata%qinfl_subsurf_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qinfl_subsurf_clms,ierr)
    if(clm_pf_idata%qsurf_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qsurf_subsurf_pfp,ierr)
    if(clm_pf_idata%qsurf_subsurf_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qsurf_subsurf_clms,ierr)
    if(clm_pf_idata%qflux_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflux_subbase_pfp,ierr)
    if(clm_pf_idata%qflux_subbase_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflux_subbase_clms,ierr)

    if(clm_pf_idata%eflux_subsurf_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subsurf_pfp,ierr)
    if(clm_pf_idata%eflux_subsurf_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subsurf_clms,ierr)
    if(clm_pf_idata%eflux_subbase_pfp /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subbase_pfp,ierr)
    if(clm_pf_idata%eflux_subbase_clms /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflux_subbase_clms,ierr)

    !-----
    if(clm_pf_idata%qflow_pfp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflow_pfp,ierr)
    if(clm_pf_idata%qflow_clms  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflow_clms,ierr)
    if(clm_pf_idata%qflowt_pfp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflowt_pfp,ierr)
    if(clm_pf_idata%qflowt_clms  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%qflowt_clms,ierr)
    if(clm_pf_idata%eflow_pfp  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflow_pfp,ierr)
    if(clm_pf_idata%eflow_clms  /= PETSC_NULL_VEC) &
       call VecDestroy(clm_pf_idata%eflow_clms,ierr)

    ! -----------------------------------------------------------------------------------------------------------

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data

#endif
