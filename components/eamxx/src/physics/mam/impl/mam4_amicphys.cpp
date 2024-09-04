#include <mam4xx/aging.hpp>
#include <mam4xx/coagulation.hpp>
#include <mam4xx/gas_chem_mechanism.hpp>
#include <mam4xx/gasaerexch.hpp>
#include <mam4xx/mam4.hpp>
#include <mam4xx/nucleation.hpp>

namespace scream::impl {

using namespace mam4;

// number of constituents in gas chemistry "work arrays"
using mam4::gas_chemistry::gas_pcnst;

// number of modes in modal configuration
constexpr int num_modes = AeroConfig::num_modes();

// number of gases
constexpr int num_gas_ids = AeroConfig::num_gas_ids();

// number of aerosol species
constexpr int num_aerosol_ids = AeroConfig::num_aerosol_ids();

//-----------------------------------------------------------------------------
//                           Indices for amicphys
//-----------------------------------------------------------------------------

// Indices of aerosol number for the arrays dimensioned gas_pcnst
KOKKOS_INLINE_FUNCTION int numptr_amode_gas_pcnst(const int mode) {
  static constexpr int numptr_amode_gas_pcnst_[num_modes] = {13, 18, 26, 30};
  return numptr_amode_gas_pcnst_[mode];
}

// Indices of aerosol mass for the arrays dimensioned gas_pcnst
KOKKOS_INLINE_FUNCTION int lmassptr_amode_gas_pcnst(const int aero_id,
                                                    const int mode) {
  static constexpr int lmassptr_amode_gas_pcnst_[num_aerosol_ids][num_modes] = {
      {6, 14, 19, 27},  {7, 15, 20, 28},  {8, 16, 21, 29}, {9, 17, 22, -6},
      {10, -6, 23, -6}, {11, -6, 24, -6}, {12, -6, 25, -6}};
  return lmassptr_amode_gas_pcnst_[aero_id][mode];
}
KOKKOS_INLINE_FUNCTION constexpr int lmapcc_val_nul() { return 0; }
KOKKOS_INLINE_FUNCTION constexpr int lmapcc_val_gas() { return 1; }
KOKKOS_INLINE_FUNCTION constexpr int lmapcc_val_aer() { return 2; }
KOKKOS_INLINE_FUNCTION constexpr int lmapcc_val_num() { return 3; }
KOKKOS_INLINE_FUNCTION int lmapcc_all(const int index) {
  static constexpr int lmapcc_all_[gas_pcnst] = {
      lmapcc_val_nul(), lmapcc_val_nul(), lmapcc_val_gas(), lmapcc_val_nul(),
      lmapcc_val_nul(), lmapcc_val_gas(), lmapcc_val_aer(), lmapcc_val_aer(),
      lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_aer(),
      lmapcc_val_aer(), lmapcc_val_num(), lmapcc_val_aer(), lmapcc_val_aer(),
      lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_num(), lmapcc_val_aer(),
      lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_aer(),
      lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_num(), lmapcc_val_aer(),
      lmapcc_val_aer(), lmapcc_val_aer(), lmapcc_val_num()};
  return lmapcc_all_[index];
}

// Where lmapcc_val_gas are defined in lmapcc_all
KOKKOS_INLINE_FUNCTION int lmap_gas(const int mode) {
  static constexpr int lmap_gas_[num_modes] = {5, 2};
  return lmap_gas_[mode];
}

// number of aerosol/gas species tendencies
KOKKOS_INLINE_FUNCTION
constexpr int nqtendbb() { return 4; }

KOKKOS_INLINE_FUNCTION constexpr int nqtendaa() { return 4; }
KOKKOS_INLINE_FUNCTION constexpr int nqqcwtendaa() { return 1; }
KOKKOS_INLINE_FUNCTION constexpr int iqtend_cond() { return 0; }
KOKKOS_INLINE_FUNCTION constexpr int iqtend_rnam() { return 1; }
KOKKOS_INLINE_FUNCTION constexpr int iqtend_nnuc() { return 2; }
KOKKOS_INLINE_FUNCTION constexpr int iqtend_coag() { return 3; }
KOKKOS_INLINE_FUNCTION constexpr int iqtend_cond_only() { return 4; }
KOKKOS_INLINE_FUNCTION constexpr int iqqcwtend_rnam() { return 0; }
KOKKOS_INLINE_FUNCTION constexpr int n_agepair() { return 1; }

// In EAMv2, subarea can take 3 values (0,1 and 2), therefore
//  length of the maxsubarea is 3
KOKKOS_INLINE_FUNCTION constexpr int maxsubarea() { return 3; }

// number of gases used in aerosol microphysics (soag and h2so4)
KOKKOS_INLINE_FUNCTION constexpr int max_gas() { return 2; }

// Index for h2so4 and nh3
constexpr int igas_h2so4 = 1;   // FIXME: This can change with modal model
constexpr int igas_nh3   = -1;  // FIXME: This can change with modal model

// leave number mix-ratios unchanged (#/kmol-air)
KOKKOS_INLINE_FUNCTION Real fcvt_num() { return 1; }

// factor for converting aerosol water mix-ratios from (kg/kg) to (mol/mol)
KOKKOS_INLINE_FUNCTION Real fcvt_wtr() {
  return haero::Constants::molec_weight_dry_air /
         haero::Constants::molec_weight_h2o;
}

// Returns index of aerosol numbers in gas_pcnst array
KOKKOS_INLINE_FUNCTION int lmap_num(const int mode) {
  return numptr_amode_gas_pcnst(mode);
}

// Returns index of aerosol numbers in gas_pcnst array
KOKKOS_INLINE_FUNCTION int lmap_numcw(const int mode) {
  return numptr_amode_gas_pcnst(mode);
}
// aerosol mapping for aerosol microphysics
// NOTE: it is different from "lmassptr_amode_gas_pcnst" as
//       amicphys adds aerosol species in a special order that is different from
//       lmassptr_amode_gas_pcnst
KOKKOS_INLINE_FUNCTION int lmap_aer(const int iaer, const int mode) {
  static constexpr int lmap_aer_[num_aerosol_ids][num_modes] = {
      {8, 15, 24, -1},  {6, 14, 21, -1},  {7, -1, 23, 27},  {9, -1, 22, 28},
      {11, 16, 20, -1}, {10, -1, 19, -1}, {12, 17, 25, 29},
  };
  return lmap_aer_[iaer][mode];
}

KOKKOS_INLINE_FUNCTION int lmap_aercw(const int iaer, const int mode) {
  return lmap_aer(iaer, mode);
}

constexpr Real mass_2_vol[num_aerosol_ids] = {0.15,
                                              6.4971751412429377e-002,
                                              0.15,
                                              7.0588235294117650e-003,
                                              3.0789473684210526e-002,
                                              5.1923076923076926e-002,
                                              156.20986883198000};

// conversion factor for aerosols
// NOTE: The following array has a special order to match amicphys
KOKKOS_INLINE_FUNCTION Real fcvt_aer(const int iaer) {
  static constexpr Real fcvt_aer_[num_aerosol_ids] = {
      8.000000000000000E-002, 1, 8.000000000000000E-002, 1, 1, 1, 1};
  return fcvt_aer_[iaer];
}

// Number of differently tagged secondary-organic aerosol species
KOKKOS_INLINE_FUNCTION constexpr int nsoa() { return 1; }

// conversion factor for gases
KOKKOS_INLINE_FUNCTION Real fcvt_gas(const int gas_id) {
  // mw to use for soa
  // BAD CONSTANTS
  constexpr Real mwuse_soa = 150;
  // molecular weight of the gas
  Real mw_gas = mam4::gas_chemistry::adv_mass[lmap_gas(gas_id)];
  // denominator
  Real denom = mw_gas;
  // special case for soa
  if(gas_id < nsoa()) denom = mwuse_soa;
  return mw_gas / denom;
}

//--------------------------------------------------------------------------------
// Utility functions
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void copy_1d_array(const int arr_len, const Real (&arr_in)[arr_len],  // in
                   Real (&arr_out)[arr_len]) {                        // out
  for(int ii = 0; ii < arr_len; ++ii) {
    arr_out[ii] = arr_in[ii];
  }
}

KOKKOS_INLINE_FUNCTION
void copy_2d_array(const int first_dimlen,                             // in
                   const int second_dimlen,                            // in
                   const Real (&arr_in)[first_dimlen][second_dimlen],  // in
                   Real (&arr_out)[first_dimlen][second_dimlen]) {     // out

  for(int ifd = 0; ifd < first_dimlen; ++ifd) {
    for(int isd = 0; isd < second_dimlen; ++isd) {
      arr_out[ifd][isd] = arr_in[ifd][isd];
    }
  }
}
template <typename DT>
KOKKOS_INLINE_FUNCTION void assign_1d_array(const int arr_len,
                                            const DT num,   // in
                                            DT *arr_out) {  // out
  for(int ii = 0; ii < arr_len; ++ii) {
    arr_out[ii] = num;
  }
}

KOKKOS_INLINE_FUNCTION
void assign_2d_array(const int first_dimlen,                          // in
                     const int second_dimlen,                         // in
                     const Real num,                                  // in
                     Real (&arr_out)[first_dimlen][second_dimlen]) {  // out
  for(int ifd = 0; ifd < first_dimlen; ++ifd) {
    for(int isd = 0; isd < second_dimlen; ++isd) {
      arr_out[ifd][isd] = num;
    }
  }
}
// copy 3d arrays
KOKKOS_INLINE_FUNCTION
void assign_3d_array(
    const int first_dimlen,                                        // in
    const int second_dimlen,                                       // in
    const int third_dimlen,                                        // in
    const Real num,                                                // in
    Real (&arr_out)[first_dimlen][second_dimlen][third_dimlen]) {  // out
  for(int ifd = 0; ifd < first_dimlen; ++ifd) {
    for(int isd = 0; isd < second_dimlen; ++isd) {
      for(int itd = 0; itd < third_dimlen; ++itd) {
        arr_out[ifd][isd][itd] = num;
      }
    }
  }
}

//--------------------------------------------------------------------------------
// Configuration settings
//--------------------------------------------------------------------------------

// MAM4 aerosol microphysics configuration data
struct AmicPhysConfig {
  // these switches activate various aerosol microphysics processes
  bool do_cond;    // condensation (a.k.a gas-aerosol exchange)
  bool do_rename;  // mode "renaming"
  bool do_newnuc;  // gas -> aerosol nucleation
  bool do_coag;    // aerosol coagulation

  // controls treatment of h2so4 condensation in mam_gasaerexch_1subarea
  //    1 = sequential   calc. of gas-chem prod then condensation loss
  //    2 = simultaneous calc. of gas-chem prod and  condensation loss
  int gaexch_h2so4_uptake_optaa;

  // controls how nucleation interprets h2so4 concentrations
  int newnuc_h2so4_conc_optaa;
};

namespace {

KOKKOS_INLINE_FUNCTION
void setup_subareas(const Real cld,                          // in
                    int &nsubarea, int &ncldy_subarea,       // out
                    int &jclea, int &jcldy,                  // out
                    bool (&iscldy_subarea)[(maxsubarea())],  // out
                    Real (&afracsub)[maxsubarea()],          // out
                    Real &fclea, Real &fcldy)                // out
{
  //--------------------------------------------------------------------------------------
  // Purpose: Determine the number of sub-areas and their fractional areas.
  //          Assign values to some bookkeeping variables.
  //--------------------------------------------------------------------------------------

  // cld: cloud fraction in the grid cell [unitless]
  // nsubarea: total # of subareas to do calculations for
  // ncldy_subarea: total # of cloudy subareas
  // jclea, jcldy: indices of the clear and cloudy subareas
  // iscldy_subarea(maxsubarea): whether a subarea is cloudy
  // afracsub(maxsubarea): area fraction of each active subarea[unitless]
  // fclea, fcldy: area fraction of clear/cloudy subarea [unitless]

  // BAD CONSTANT
  //  Cloud chemistry is only active when cld(i,kk) >= 1.0e-5
  //  It may be that the macrophysics has a higher threshold than this
  constexpr Real fcld_locutoff = 1.0e-5;

  // BAD CONSTANT
  //  Grid cells with cloud fraction larger than this cutoff is considered to be
  //  overcast
  constexpr Real fcld_hicutoff = 0.999;

  // if cloud fraction ~= 0, the grid-cell has a single clear  sub-area
  // (nsubarea = 1) if cloud fraction ~= 1, the grid-cell has a single cloudy
  // sub-area (nsubarea = 1) otherwise, the grid-cell has a clear and a cloudy
  // sub-area (nsubarea = 2)

  if(cld < fcld_locutoff) {
    fcldy         = 0;
    nsubarea      = 1;
    ncldy_subarea = 0;
    jclea         = 1;
    jcldy         = 0;
  } else if(cld > fcld_hicutoff) {
    fcldy         = 1;
    nsubarea      = 1;
    ncldy_subarea = 1;
    jclea         = 0;
    jcldy         = 1;
  } else {
    fcldy         = cld;
    nsubarea      = 2;
    ncldy_subarea = 1;
    jclea         = 1;
    jcldy         = 2;
  }

  fclea = 1.0 - fcldy;

  // Set up a logical array to indicate whether the subareas are clear or cloudy
  // and init area fraction array
  for(int jsub = 0; jsub < maxsubarea(); ++jsub) {
    iscldy_subarea[jsub] = false;
    afracsub[jsub]       = 0;
  }

  // jcldy>0 can be 1 or 2, so iscldy_subarea(1) or iscldy_subarea(2) is true
  if(jcldy > 0) iscldy_subarea[jcldy] = true;
  // Save the area fractions to an array
  // jclea can only be 1 if jclea > 0, so afracsub (1) is set to fclea
  if(jclea > 0) afracsub[jclea] = fclea;
  // jcldy can be 1 or 2, so afracsub(1) or afracsub(2) is set to fcldy
  if(jcldy > 0) afracsub[jcldy] = fcldy;

}  // setup_subareas

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void set_subarea_rh(const int &ncldy_subarea, const int &jclea,  // in
                    const int &jcldy,                            // in
                    const Real (&afracsub)[maxsubarea()],        // in
                    const Real &relhumgcm,                       // in
                    Real (&relhumsub)[maxsubarea()])             // out
{
  //----------------------------------------------------------------------------
  // Purpose: Set relative humidity in subareas.
  //----------------------------------------------------------------------------

  // ncldy_subarea         :# of cloudy subareas
  // jclea, jcldy          :indices of clear and cloudy subareas
  // afracsub(maxsubarea)  :area fraction in subareas [unitless]
  // relhumgcm             :grid cell mean relative humidity [unitless]
  // relhumsub(maxsubarea): relative humidity in subareas [unitless]

  if(ncldy_subarea <= 0) {
    // Entire grid cell is cloud-free. RH in subarea = grid cell mean.
    // This is clear cell, rehumsub(0),rehumsub(1) and rehumsub(3) are relhumgcm
    for(int jsub = 0; jsub < maxsubarea(); ++jsub) relhumsub[jsub] = relhumgcm;
  } else {
    // Grid cell has a cloudy subarea. Set RH in that part to 1.0.
    // jcldy can be 1 or 2 here.
    // If jcldy is 1: relhumsub[1] is 1.0 (fully cloudy cell)
    // if jcldy is 2: relhumsub[2] is 1.0. In this case jclea is >0,
    //               so relhumsub[1] is set in if condition below
    relhumsub[jcldy] = 1;

    // If the grid cell also has a clear portion, back out the subarea RH from
    // the grid-cell mean RH and cloud fraction.
    if(jclea > 0) {
      // jclea is > 0 only for partly cloudy cell. In this case
      // jclea is 1, so relhumsub[1] is set here.
      Real relhum_tmp  = (relhumgcm - afracsub[jcldy]) / afracsub[jclea];
      relhumsub[jclea] = mam4::utils::min_max_bound(0, 1, relhum_tmp);
    }
  }
}  // set_subarea_rh

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void compute_qsub_from_gcm_and_qsub_of_other_subarea(
    const bool (&lcompute)[gas_pcnst], const Real &f_a,  // in
    const Real &f_b,                                     // in
    const Real (&qgcm)[gas_pcnst],                       // in
    const int &jclea, const int &jcldy,                  // in
    Real (&qsub_a)[gas_pcnst][maxsubarea()],             // inout
    Real (&qsub_b)[gas_pcnst][maxsubarea()])             // inout
{
  //-----------------------------------------------------------------------------------------
  // Purpose: Calculate the value of qsub_b assuming qgcm is a weighted average
  // defined as
  //          qgcm = f_a*qsub_a + f_b*qsub_b.
  //-----------------------------------------------------------------------------------------

  // f_a, f_b      // area fractions [unitless] of subareas
  // qgcm(ncnst)   // grid cell mean (known)
  // qsub_a(ncnst) // value in subarea A (known, but might get adjusted)
  // qsub_b(ncnst) // value in subarea B (to be calculated here)

  // Here we populate qsub for subarea index 2 (i.e. jcldy is 2 here)
  //  and adjust subarea index 1(i.e., jclea is 1 here) if needed.
  for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
    if(lcompute[icnst]) {
      // Calculate qsub_b
      EKAT_KERNEL_ASSERT_MSG(
          f_b != 0,
          "Error! compute_qsub_from_gcm_and_qsub_of_other_subarea - f_b is "
          "zero\n");
      qsub_b[icnst][jcldy] = (qgcm[icnst] - f_a * qsub_a[icnst][jclea]) / f_b;
      // Check that this does not produce a negative value.
      // If so, set qsub_b to zero and adjust the value of qsub_a.
      if(qsub_b[icnst][jcldy] < 0) {
        qsub_b[icnst][jcldy] = 0;
        EKAT_KERNEL_ASSERT_MSG(
            f_a != 0,
            "Error! compute_qsub_from_gcm_and_qsub_of_other_subarea - f_a is "
            "zero\n");
        qsub_a[icnst][jclea] = qgcm[icnst] / f_a;
      }
    }
  }
}  // compute_qsub_from_gcm_and_qsub_of_other_subarea

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void set_subarea_qnumb_for_cldbrn_aerosols(
    const int &jclea, const int &jcldy, const Real &fcldy,  // in
    const Real (&qqcwgcm)[gas_pcnst],                       // in
    Real (&qqcwsub)[gas_pcnst][maxsubarea()])               // inout
{
  //-----------------------------------------------------------------------------------------
  // Purpose: Set the number mixing ratios of cloud-borne aerosols in subareas:
  //          - zero in clear air;
  //          - grid-cell-mean divided by cloud-fraction in the cloudy subarea.
  //          This is done for all lognormal modes.
  //-----------------------------------------------------------------------------------------

  // jclea, jcldy              : indices of subareas
  // fcldy                     : area fraction [unitless] of the cloudy subarea
  // qqcwgcm(ncnst)            : grid cell mean (unit does not matter for this
  //                             subr.)
  // qqcwsub(ncnst,maxsubarea) : values in subareas (unit does not matter
  //                             for this subr.)

  //----------------------------------------------------------------
  // Here jclea ==1 and jcldy==2
  for(int imode = 0; imode < num_modes; ++imode) {
    const int icnst       = numptr_amode_gas_pcnst(imode);
    qqcwsub[icnst][jclea] = 0;
    EKAT_KERNEL_ASSERT_MSG(
        fcldy != 0,
        "Error! set_subarea_qnumb_for_cldbrn_aerosols - fcldy is "
        "zero\n");
    qqcwsub[icnst][jcldy] = qqcwgcm[icnst] / fcldy;
    //----------------------------------------------------------------
  }

}  // set_subarea_qnumb_for_cldbrn_aerosols

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void set_subarea_qmass_for_cldbrn_aerosols(
    const int &jclea, const int &jcldy,        // in
    const Real &fcldy,                         // in
    const Real (&qqcwgcm)[gas_pcnst],          // in
    Real (&qqcwsub)[gas_pcnst][maxsubarea()])  // inout
{
  //-----------------------------------------------------------------------------------------
  // Purpose: Set the mass mixing ratios of cloud-borne aerosols in subareas:
  //          - zero in clear air;
  //          - grid-cell-mean/cloud-fraction in the cloudy subarea.
  //          This is done for all lognormal modes and all chemical species.
  //-----------------------------------------------------------------------------------------
  // jclea, jcldy              : subarea indices fcldy : area
  //                             fraction [unitless] of the cloudy subarea
  // qqcwgcm(ncnst)            : grid cell mean (unit does not matter for this
  //                             subr.)
  // qqcwsub(ncnst,maxsubarea) : values in subareas (unit does not matter for
  //                             this subr.)

  //----------------------------------------------------------------
  // Here jclea ==1 and jcldy==2

  // loop thru all modes
  for(int imode = 0; imode < num_modes; ++imode) {
    // loop thru all species in a mode
    for(int ispec = 0; ispec < mam4::num_species_mode(imode); ++ispec) {
      const int icnst = lmassptr_amode_gas_pcnst(ispec, imode);

      qqcwsub[icnst][jclea] = 0;
      EKAT_KERNEL_ASSERT_MSG(
          fcldy != 0,
          "Error! set_subarea_qmass_for_cldbrn_aerosols - fcldy is "
          "zero\n");
      qqcwsub[icnst][jcldy] = qqcwgcm[icnst] / fcldy;
    }  // ispec - species loop
  }    // imode - mode loop
       //----------------------------------------------------------------
}  // set_subarea_qmass_for_cldbrn_aerosols

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void get_partition_factors(const Real &qgcm_intrst,               // in
                           const Real &qgcm_cldbrn,               // in
                           const Real &fcldy, const Real &fclea,  // in
                           Real &factor_clea, Real &factor_cldy)  // out
{
  //------------------------------------------------------------------------------------
  // Purpose: Calculate the partitioning factors for distributing interstitial
  // aerosol
  //          mixing ratios to cloudy and clear subareas in a grid box.
  //          The partitioning factors depend on the grid cell mean mixing
  //          ratios of both interstitial and cloud-borne aerosols.
  //------------------------------------------------------------------------------------

  // qgcm_intrst  : grid cell mean interstitial aerosol mixing ratio
  // qgcm_cldbrn  : grid cell mean cloud-borne aerosol mixing ratio
  //
  // fcldy        : cloudy fraction of the grid cell [unitless]
  // fclea        : clear  fraction of the grid cell [unitless]
  //
  // factor_clea  : partitioning factor for clear  subarea
  // factor_cldy  : partitioning factor for cloudy subarea

  // Calculate subarea-mean mixing ratios

  EKAT_KERNEL_ASSERT_MSG(fcldy != 0,
                         "Error! get_partition_factors - fcldy is "
                         "zero\n");
  // cloud-borne,  cloudy subarea
  const Real tmp_q_cldbrn_cldy = qgcm_cldbrn / fcldy;

  // interstitial, cloudy subarea
  const Real tmp_q_intrst_cldy =
      haero::max(0, ((qgcm_intrst + qgcm_cldbrn) - tmp_q_cldbrn_cldy));

  EKAT_KERNEL_ASSERT_MSG(fclea != 0,
                         "Error! get_partition_factors - fclea is "
                         "zero\n");
  // interstitial, clear  subarea
  const Real tmp_q_intrst_clea =
      (qgcm_intrst - fcldy * tmp_q_intrst_cldy) / fclea;

  // Calculate the corresponding paritioning factors for interstitial
  // aerosols using the above-derived subarea-mean mixing ratios plus the
  // constraint that the cloud fraction weighted average of subarea mean
  // need to match grid box mean. Note that this subroutine is designed for
  // partially cloudy grid cells, hence both fclea and fcldy are assumed to
  // be nonzero.

  constexpr Real eps = 1.e-35;  // BAD CONSTANT
  Real clea2gcm_ratio =
      haero::max(eps, tmp_q_intrst_clea * fclea) / haero::max(eps, qgcm_intrst);
  clea2gcm_ratio = haero::max(0, haero::min(1, clea2gcm_ratio));

  factor_clea = clea2gcm_ratio / fclea;
  factor_cldy = (1 - clea2gcm_ratio) / fcldy;
}  // get_partition_factors

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void set_subarea_qnumb_for_intrst_aerosols(
    const int &jclea, const int &jcldy, const Real &fclea,  // in
    const Real &fcldy, const Real (&qgcm)[gas_pcnst],       // in
    const Real (&qqcwgcm)[gas_pcnst],                       // in
    const Real (&qgcmx)[gas_pcnst],                         // in
    Real (&qsubx)[gas_pcnst][maxsubarea()])                 // inout
{
  //-----------------------------------------------------------------------------------------
  // Purpose: Set the number mixing ratios of interstitial aerosols in subareas.
  //          Interstitial aerosols can exist in both cloudy and clear subareas,
  //          so a grid cell mean needs to be partitioned. Different lognormal
  //          modes are partitioned differently based on the mode-specific
  //          number mixing ratios.
  //-----------------------------------------------------------------------------------------

  // jclea, jcldy  : subarea indices
  // fclea, fcldy  : area fraction [unitless] of the clear and cloudy subareas
  // qgcm   (ncnst): grid cell mean, interstitial constituents (unit does not
  //                 matter)
  // qqcwgcm(ncnst): grid cell mean, cloud-borne  constituents (unit
  //                 does not matter)

  // qgcmx  (ncnst): grid cell mean, interstitial constituents (unit does not
  //                 matter)
  // qsubx(ncnst,maxsubarea): subarea mixing ratios of interst. constituents
  //                          (unit does not matter as long as they are
  //                          consistent with the grid cell mean values)

  // Note: qgcm and qqcwgcm are used for calculating the patitioning factors.
  // qgcmx is the actual grid cell mean that is partitioned into qsubx.

  for(int imode = 0; imode < num_modes; ++imode) {
    // calculate partitioning factors

    // grid cell mean of interstitial aerosol mixing ratio of a single mode
    const Real qgcm_intrst = qgcm[numptr_amode_gas_pcnst(imode)];

    // grid cell mean of cloud-borne  aerosol mixing ratio of a single mode
    const Real qgcm_cldbrn = qqcwgcm[numptr_amode_gas_pcnst(imode)];

    Real factor_clea;  // partitioning factor for clear  subarea [unitless]
    Real factor_cldy;  // partitioning factor for cloudy subarea [unitless]
    get_partition_factors(qgcm_intrst, qgcm_cldbrn, fcldy, fclea,  // in
                          factor_clea, factor_cldy);               // out

    // apply partitioning factors
    const int icnst = numptr_amode_gas_pcnst(imode);

    qsubx[icnst][jclea] = qgcmx[icnst] * factor_clea;
    qsubx[icnst][jcldy] = qgcmx[icnst] * factor_cldy;
  }  // imode

}  // set_subarea_qnumb_for_intrst_aerosols

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void set_subarea_qmass_for_intrst_aerosols(
    const int &jclea, const int &jcldy, const Real &fclea,  // in
    const Real &fcldy, const Real (&qgcm)[gas_pcnst],       // in
    const Real (&qqcwgcm)[gas_pcnst],                       // in
    const Real (&qgcmx)[gas_pcnst],                         // in
    Real (&qsubx)[gas_pcnst][maxsubarea()])                 // inout
{
  //-----------------------------------------------------------------------------------------
  // Purpose: Set the mass mixing ratios of interstitial aerosols in subareas.
  //          Interstitial aerosols can exist in both cloudy and clear subareas,
  //          so a grid cell mean needs to be partitioned. Different lognormal
  //          modes are partitioned differently based on the mode-specific
  //          mixing ratios. All species in the same mode are partitioned the
  //          same way, consistent with the internal mixing assumption used in
  //          MAM.
  //-----------------------------------------------------------------------------------------

  // jclea, jcldy   : subarea indices
  // fclea, fcldy   : area fraction [unitless] of the clear and cloudy subareas
  // qgcm   (ncnst) : grid cell mean, interstitial constituents (unit does not
  //                  matter)
  // qqcwgcm(ncnst) : grid cell mean, cloud-borne  constituents (unit
  //                  does not matter)

  // qgcmx  (ncnst) : grid cell mean, interstitial constituents (unit does not
  //                  matter)
  // qsubx(ncnst,maxsubarea): subarea mixing ratios of interst.
  //                          constituents(unit does not matter as long as they
  //                          are consistent with the grid cell mean values)

  // Note: qgcm and qqcwgcm are used for calculating the patitioning factors.
  // qgcmx is the actual grid cell mean that is partitioned into qsubx.

  for(int imode = 0; imode < num_modes; ++imode) {
    // calculcate partitioning factors

    // grid cell mean of interstitial aerosol mixing ratio of a single mode
    Real qgcm_intrst = 0;

    // grid cell mean of cloud-borne  aerosol mixing ratio of a single mode
    Real qgcm_cldbrn = 0;

    // loop thru all species in a mode
    for(int ispec = 0; ispec < mam4::num_species_mode(imode); ++ispec) {
      qgcm_intrst = qgcm_intrst + qgcm[lmassptr_amode_gas_pcnst(ispec, imode)];
      qgcm_cldbrn =
          qgcm_cldbrn + qqcwgcm[lmassptr_amode_gas_pcnst(ispec, imode)];
    }

    Real factor_clea;  // partitioning factor for clear  subarea [unitless]
    Real factor_cldy;  // partitioning factor for cloudy subarea [unitless]
    get_partition_factors(qgcm_intrst, qgcm_cldbrn, fcldy, fclea,  // in
                          factor_clea, factor_cldy);               // out

    // apply partitioning factors
    // Here jclea==1 and jcldy==2
    for(int ispec = 0; ispec < mam4::num_species_mode(imode); ++ispec) {
      const int icnst     = lmassptr_amode_gas_pcnst(ispec, imode);
      qsubx[icnst][jclea] = qgcmx[icnst] * factor_clea;
      qsubx[icnst][jcldy] = qgcmx[icnst] * factor_cldy;
    }  // ispec
  }    // imode

}  // set_subarea_qmass_for_intrst_aerosols

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void set_subarea_gases_and_aerosols(
    const int &nsubarea, const int &jclea,                           // in
    const int &jcldy,                                                // in
    const Real &fclea, const Real &fcldy,                            // in
    const Real (&qgcm1)[gas_pcnst], const Real (&qgcm2)[gas_pcnst],  // in
    const Real (&qqcwgcm2)[gas_pcnst],                               // in
    const Real (&qgcm3)[gas_pcnst],                                  // in
    const Real (&qqcwgcm3)[gas_pcnst],                               // in
    Real (&qsub1)[gas_pcnst][maxsubarea()],                          // out
    Real (&qsub2)[gas_pcnst][maxsubarea()],                          // out
    Real (&qqcwsub2)[gas_pcnst][maxsubarea()],                       // out
    Real (&qsub3)[gas_pcnst][maxsubarea()],                          // out
    Real (&qqcwsub3)[gas_pcnst][maxsubarea()])                       // out
{
  //------------------------------------------------------------------------------------------------
  // Purpose: Partition grid cell mean mixing ratios to clear/cloudy subareas.
  //------------------------------------------------------------------------------------------------
  // nsubarea: # of active subareas in the current grid cell
  // jclea, jcldy: indices of the clear and cloudy subareas
  // fclea, fcldy: area fractions of the clear and cloudy subareas [unitless]

  // The next set of argument variables are tracer mixing ratios.
  //  - The units are different for gases, aerosol number, and aerosol mass.
  //    The exact units do not matter for this subroutine, as long as the
  //    grid cell mean values ("gcm") and the corresponding subarea values
  //    ("sub") have the same units.
  //  - qq* and qqcw* are correspond to the interstitial and cloud-borne
  //    species, respectively
  //  - The numbers 1-3 correspond to different locations in the host model's
  //    time integration loop.

  // Grid cell mean mixing ratios
  // qgcm1(ncnst), qgcm2(ncnst), qqcwgcm2(ncnst), qgcm3(ncnst),
  // qqcwgcm3(ncnst)

  // Subarea mixing ratios
  // qsub1(ncnst,maxsubarea), qsub2(ncnst,maxsubarea),
  // qsub3(ncnst,maxsubarea), qqcwsub2(ncnst,maxsubarea)
  // qqcwsub3(ncnst,maxsubarea)
  //----

  //------------------------------------------------------------------------------------
  // Initialize mixing ratios in subareas before the aerosol microphysics
  // calculations
  //------------------------------------------------------------------------------------
  // FIXME:Should we set jsub==0 a special value (like NaNs) so that it is never
  // used??
  for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
    for(int jsub = 0; jsub < maxsubarea(); ++jsub) {
      // Gases and interstitial aerosols
      qsub1[icnst][jsub] = 0;
      qsub2[icnst][jsub] = 0;
      qsub3[icnst][jsub] = 0;

      // Cloud-borne aerosols
      qqcwsub2[icnst][jsub] = 0;
      qqcwsub3[icnst][jsub] = 0;
    }
  }
  //---------------------------------------------------------------------------------------------------
  // Determine which category the current grid cell belongs to: partly cloudy,
  // all cloudy, or all clear
  //---------------------------------------------------------------------------------------------------
  const bool grid_cell_has_only_clea_area =
      ((jclea == 1) && (jcldy == 0) && (nsubarea == 1));
  const bool grid_cell_has_only_cldy_area =
      ((jclea == 0) && (jcldy == 1) && (nsubarea == 1));
  const bool gird_cell_is_partly_cldy =
      (jclea > 0) && (jcldy > 0) && (jclea + jcldy == 3) && (nsubarea == 2);

  // Sanity check
  if((!grid_cell_has_only_clea_area) && (!grid_cell_has_only_cldy_area) &&
     (!gird_cell_is_partly_cldy)) {
    EKAT_KERNEL_ASSERT_MSG(true,
                           "Error! modal_aero_amicphys - bad jclea, jcldy, "
                           "nsubarea, jclea, jcldy, nsubarea\n");
  }

  //*************************************************************************************************
  // Category I: grid cell is either all clear or all cloudy. Copy the grid
  // cell mean values.
  //*************************************************************************************************
  if(grid_cell_has_only_clea_area || grid_cell_has_only_cldy_area) {
    // For fully clear and cloudy cells, we populate only 1st index of subarea
    // for all output vars
    //  Makes sense as there is only 1 subarea for these cases.
    // FIXME: Should we fill in NaNs for the 0th and 2nd index??
    constexpr int jsub = 1;
    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      // copy all gases and aerosols
      if(lmapcc_all(icnst) > 0) {
        qsub1[icnst][jsub] = qgcm1[icnst];
        qsub2[icnst][jsub] = qgcm2[icnst];
        qsub3[icnst][jsub] = qgcm3[icnst];

        qqcwsub2[icnst][jsub] = qqcwgcm2[icnst];
        qqcwsub3[icnst][jsub] = qqcwgcm3[icnst];
      }
    }
  }  // if only clear or only cloudy
  //*************************************************************************************************
  // Category II: partly cloudy grid cell. Tracer mixing ratios are generally
  // assumed different in clear and cloudy subareas.  This is primarily
  // because the interstitial aerosol mixing ratios are assumed to be lower
  // in the cloudy sub-area than in the clear sub-area, as much of the
  // aerosol is activated in the cloudy sub-area.
  //*************************************************************************************************
  else if(gird_cell_is_partly_cldy) {
    //===================================
    // Set gas mixing ratios in subareas
    //===================================
    //------------------------------------------------------------------------------------------
    // Before gas chemistry, gas mixing ratios are assumed to be the same in
    // all subareas, i.e., they all equal the grid cell mean.
    //------------------------------------------------------------------------------------------

    // NOTE: In this "else if" case jclea == 1 and jcldy == 2

    bool cnst_is_gas[gas_pcnst] = {};
    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      cnst_is_gas[icnst] = (lmapcc_all(icnst) == lmapcc_val_gas());
    }

    EKAT_KERNEL_ASSERT_MSG(nsubarea < maxsubarea(),
                           "Error! set_subarea_gases_and_aerosols: "
                           "nsubarea should be < maxsubarea() \n");
    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      if(cnst_is_gas[icnst]) {
        // For gases, assume both 1 and 2 subareas have grid mean values
        for(int jsub = 1; jsub <= nsubarea; ++jsub) {
          qsub1[icnst][jsub] = qgcm1[icnst];
        }
      }
    }
    // qsub1 is fully populated for gasses
    //------------------------------------------------------------------------------------------
    // After gas chemistry, still assume gas mixing ratios are the same in all
    // subareas.
    //------------------------------------------------------------------------------------------

    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      if(cnst_is_gas[icnst]) {
        // For gases, assume both 1 and 2 subareas have grid mean values
        for(int jsub = 1; jsub <= nsubarea; ++jsub) {
          qsub2[icnst][jsub] = qgcm2[icnst];
        }
      }
    }
    //  qsub2 is fully populated for gasses
    //----------------------------------------------------------------------------------------
    //   After cloud chemistry, gas and aerosol mass mixing ratios in the clear
    //   subarea are assumed to be the same as their values before cloud
    //   chemistry (because by definition, cloud chemistry did not happen in
    //   clear sky), while the mixing ratios in the cloudy subarea likely have
    //   changed.
    //----------------------------------------------------------------------------------------
    //   Gases in the clear subarea remain the same as their values before cloud
    //   chemistry.
    //  Here we populate qsub3 for index 1 only as jclea is 1.
    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      if(cnst_is_gas[icnst]) {
        qsub3[icnst][jclea] = qsub2[icnst][jclea];
      }
    }

    // Calculate the gas mixing ratios in the cloudy subarea using the
    // grid-cell mean, cloud fraction and the clear-sky values
    // Here we populate qsub3 for index 2 (jcldy) and adjust index 1 (jclea) if
    // needed.
    compute_qsub_from_gcm_and_qsub_of_other_subarea(cnst_is_gas, fclea,   // in
                                                    fcldy, qgcm3, jclea,  // in
                                                    jcldy,                // in
                                                    qsub3, qsub3);  // inout
    // qsub3[2][2]);
    //  qsub3 is fully populated for gasses
    //=========================================================================
    //  Set AEROSOL mixing ratios in subareas.
    //  Only need to do this for points 2 and 3 in the time integraion loop,
    //  i.e., the before-cloud-chem and after-cloud-chem states.
    //=========================================================================
    //  Cloud-borne aerosols. (They are straightforward to partition,
    //  as they only exist in the cloudy subarea.)
    //----------------------------------------------------------------------------------------
    //  Partition mass and number before cloud chemistry
    //  NOTE that in this case jclea is 1 and jcldy is 2
    //  Following 2 calls set qqcwsub2(:,1)=0 and qqcwsub2(:,2) to a computed
    //  value
    set_subarea_qnumb_for_cldbrn_aerosols(jclea, jcldy, fcldy,
                                          qqcwgcm2,   // in
                                          qqcwsub2);  // inout

    set_subarea_qmass_for_cldbrn_aerosols(jclea, jcldy, fcldy,
                                          qqcwgcm2,   // in
                                          qqcwsub2);  // inout
    //  Partition mass and number before cloud chemistry
    // Following 2 calls set qqcwsub3(:,1)=0 and qqcwsub3(:,2) to a computed
    // value
    set_subarea_qnumb_for_cldbrn_aerosols(jclea, jcldy, fcldy,
                                          qqcwgcm3,   // in
                                          qqcwsub3);  // inout
    set_subarea_qmass_for_cldbrn_aerosols(jclea, jcldy, fcldy,
                                          qqcwgcm3,   // in
                                          qqcwsub3);  // inout

    //----------------------------------------------------------------------------------------
    // Interstitial aerosols. (They can exist in both cloudy and clear
    // subareas, and hence need to be partitioned.)
    //----------------------------------------------------------------------------------------
    // Partition mass and number before cloud chemistry
    // Following 2 calls set qsub2(:,1) = 0 and qsub2(:,2) to a computed value
    set_subarea_qnumb_for_intrst_aerosols(jclea, jcldy, fclea, fcldy,  // in
                                          qgcm2, qqcwgcm2, qgcm2,      // in
                                          qsub2);                      // inout

    set_subarea_qmass_for_intrst_aerosols(jclea, jcldy, fclea, fcldy,  // in
                                          qgcm2, qqcwgcm2, qgcm2,      // in
                                          qsub2);                      // inout

    // Partition mass and number before cloud chemistry
    // Following 2 calls set qsub3(:,1) = 0 and qsub3(:,2) to a computed value
    set_subarea_qnumb_for_intrst_aerosols(jclea, jcldy, fclea, fcldy,  // in
                                          qgcm2, qqcwgcm2, qgcm3,      // in
                                          qsub3);                      // inout

    set_subarea_qmass_for_intrst_aerosols(jclea, jcldy, fclea, fcldy,  // in
                                          qgcm2, qqcwgcm2, qgcm3,      // in
                                          qsub3);                      // inout

  }  // different categories
}  // set_subarea_gases_and_aerosols

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void mam_newnuc_1subarea(
    const int igas_h2so4, const int gaexch_h2so4_uptake_optaa,   // in
    const int newnuc_h2so4_conc_optaa, const int jsub,           // in
    const Real deltat, const Real temp, const Real pmid,         // in
    const Real aircon, const Real zmid, const Real pblh,         // in
    const Real relhum, const Real uptkrate_h2so4,                // in
    const Real del_h2so4_gasprod, const Real del_h2so4_aeruptk,  // in
    Real qgas_cur[max_gas()], const Real qgas_avg[max_gas()],    // out
    Real qnum_cur[AeroConfig::num_modes()],                      // out
    Real qaer_cur[AeroConfig::num_aerosol_ids()]
                 [AeroConfig::num_modes()],  // out
    Real dnclusterdt) {                      // out
  // FIXME: This function was not refactored or cleaned in fortran
  //  we must clean it and remove unused codes and fix variable names

  // nstep: model time-step number
  // jsub: sub-area index
  // deltat: model timestep (s)
  // temp: temperature (K)
  // pmid: pressure at model levels(Pa)
  // aircon: air molar concentration (kmol/m3)
  // zmid: midpoint height above surface (m)
  // pblh  :pbl height (m)
  // relhum:relative humidity (0-1)
  // uptkrate_h2so4
  // del_h2so4_gasprod
  // del_h2so4_aeruptk
  // dnclusterdt: cluster nucleation rate (#/m3/s)
  // qgas_cur(max_gas)
  // qgas_avg(max_gas)
  // qnum_cur(max_mode)
  // qaer_cur(1:max_aer,1:max_mode)
  // qwtr_cur(1:max_mode)

  // DESCRIPTION:
  //   computes changes due to aerosol nucleation (new particle formation)
  //       treats both nucleation and subsequent growth of new particles
  //          to aitken mode size
  //   uses the following parameterizations
  //       vehkamaki et al. (2002) parameterization for binary
  //           homogeneous nucleation (h2so4-h2o) plus
  //       kerminen and kulmala (2002) parameterization for
  //           new particle loss during growth to aitken size
  //
  // REVISION HISTORY:
  //   R.Easter 2007.09.14:  Adapted from MIRAGE2 code and CMAQ V4.6 code
  //

  constexpr int newnuc_method_flagaa = 11;
  //  1=merikanto et al (2007) ternary   2=vehkamaki et al (2002) binary
  // 11=merikanto ternary + first-order boundary layer
  // 12=merikanto ternary + second-order boundary layer

  // begin
  dnclusterdt = 0;

  // qh2so4_cur = current qh2so4, after aeruptk
  // qh2so4_avg = average qh2so4 over time-step
  // BAD CONSTANTS
  constexpr Real qh2so4_cutoff = 4.0e-16;
  Real qh2so4_cur              = qgas_cur[igas_h2so4];
  Real qh2so4_avg, tmp_uptkrate, tmpa;

  // use qh2so4_avg and first-order loss rate calculated in
  // mam_gasaerexch_1subarea
  qh2so4_avg   = qgas_avg[igas_h2so4];
  tmp_uptkrate = uptkrate_h2so4;

  if(qh2so4_avg <= qh2so4_cutoff) return;

  static constexpr int igas_nh3 = -999888777;  // Same as mam_refactor
  Real qnh3_cur                 = 0;

  //   dry-diameter limits for "grown" new particles
  constexpr int nait = static_cast<int>(ModeIndex::Aitken);
  Real dplom_mode    = haero::exp(0.67 * haero::log(modes(nait).min_diameter) +
                                  0.33 * haero::log(modes(nait).nom_diameter));
  Real dphim_mode    = modes(nait).max_diameter;

  //   mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
  //                (assuming same dry density for so4 & nh4)
  //      mass1p_aitlo - dp = dplom_mode
  //      mass1p_aithi - dp = dphim_mode
  constexpr Real dens_so4a_host = 1770;
  tmpa                          = dens_so4a_host * haero::Constants::pi / 6.0;
  Real mass1p_aitlo             = tmpa * (haero::pow(dplom_mode, 3.0));
  Real mass1p_aithi             = tmpa * (haero::pow(dphim_mode, 3.0));

  //   limit RH to between 0.1% and 99%
  Real relhumnn = haero::max(0.01, haero::min(0.99, relhum));

  // BAD CONSTANTS (These should come from chemistry mechanism
  // but it is fixed here for stay BFB)
  constexpr Real mw_so4a_host = 115;
  constexpr Real mwnh4        = 18;
  constexpr Real mwso4        = 96;
  // BAD CONSTANTS (These should come from Haero constants
  // but it is fixed here for stay BFB)
  constexpr Real rgas     = 8.31446759100000;
  constexpr Real avogadro = 6.022140000000000E+023;

  int itmp;
  Real qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a;

  //   call ... routine to get nucleation rates
  // FIXME: I GOT ALL ZEROS...THIS IS NOT VALIDATED YET!!!!
  mam4::nucleation::mer07_veh02_nuc_mosaic_1box(
      newnuc_method_flagaa, deltat, temp, relhumnn, pmid, zmid, pblh,   // in
      qh2so4_cur, qh2so4_avg, qnh3_cur, tmp_uptkrate, mw_so4a_host, 1,  // in
      dplom_mode, dphim_mode, rgas, avogadro, mwnh4, mwso4,             // in
      haero::Constants::pi,                                             // in
      itmp, qnuma_del, qso4a_del, qnh4a_del, qh2so4_del,                // out
      qnh3_del, dens_nh4so4a, dnclusterdt);                             // out

  //   convert qnuma_del from (#/mol-air) to (#/kmol-air)
  qnuma_del = qnuma_del * 1.0e3;

  //   number nuc rate (#/kmol-air/s) from number nuc amt
  Real dndt_ait = qnuma_del / deltat;

  //   fraction of mass nuc going to so4
  tmpa           = qso4a_del * mw_so4a_host;
  Real tmpb      = tmpa;
  Real tmp_frso4 = 1.0;

  //   mass nuc rate (kg/kmol-air/s) from mass nuc amts
  EKAT_KERNEL_ASSERT_MSG(deltat != 0,
                         "Error! mam_newnuc_1subarea: "
                         " deltat should not be equal to 0\n");
  Real dmdt_ait = haero::max(0.0, (tmpb / deltat));

  Real dndt_aitsv2 = 0.0;
  Real dmdt_aitsv2 = 0.0;
  Real dndt_aitsv3 = 0.0;
  Real dmdt_aitsv3 = 0.0;
  // BAD CONSTANTS
  if(dndt_ait < 1.0e2) {
    //   ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
    dndt_ait = 0.0;
    dmdt_ait = 0.0;
  } else {
    dndt_aitsv2 = dndt_ait;
    dmdt_aitsv2 = dmdt_ait;

    //   mirage2 code checked for complete h2so4 depletion here,
    //   but this is now done in mer07_veh02_nuc_mosaic_1box
    EKAT_KERNEL_ASSERT_MSG(dndt_ait != 0,
                           "Error! mam_newnuc_1subarea: "
                           " dndt_ait should not be equal to 0\n");
    Real mass1p = dmdt_ait / dndt_ait;
    dndt_aitsv3 = dndt_ait;
    dmdt_aitsv3 = dmdt_ait;

    EKAT_KERNEL_ASSERT_MSG(mass1p_aitlo != 0,
                           "Error! mam_newnuc_1subarea: "
                           " mass1p_aitlo should not be equal to 0\n");
    //   apply particle size constraints
    if(mass1p < mass1p_aitlo) {
      //   reduce dndt to increase new particle size
      dndt_ait = dmdt_ait / mass1p_aitlo;
    } else if(mass1p > mass1p_aithi) {
      //   reduce dmdt to decrease new particle size
      dmdt_ait = dndt_ait * mass1p_aithi;
    }
  }

  // *** apply adjustment factor to avoid unrealistically high
  //     aitken number concentrations in mid and upper troposphere
  constexpr Real newnuc_adjust_factor_dnaitdt = 1;
  dndt_ait = dndt_ait * newnuc_adjust_factor_dnaitdt;
  dmdt_ait = dmdt_ait * newnuc_adjust_factor_dnaitdt;

  Real tmp_q_del = dndt_ait * deltat;
  qnum_cur[nait] = qnum_cur[nait] + tmp_q_del;

  //   dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)

  constexpr Real mw_nh4a_host = mw_so4a_host;
  EKAT_KERNEL_ASSERT_MSG(mw_so4a_host != 0,
                         "Error! mam_newnuc_1subarea: "
                         " mw_so4a_host should not be equal to 0\n");
  Real dso4dt_ait = dmdt_ait * tmp_frso4 / mw_so4a_host;
  EKAT_KERNEL_ASSERT_MSG(mw_nh4a_host != 0,
                         "Error! mam_newnuc_1subarea: "
                         " mw_nh4a_host should not be equal to 0\n");
  Real dnh4dt_ait        = dmdt_ait * (1.0 - tmp_frso4) / mw_nh4a_host;
  constexpr int iaer_so4 = 1;
  if(dso4dt_ait > 0.0) {
    tmp_q_del                = dso4dt_ait * deltat;
    qaer_cur[iaer_so4][nait] = qaer_cur[iaer_so4][nait] + tmp_q_del;
    tmp_q_del                = haero::min(tmp_q_del, qgas_cur[igas_h2so4]);
    qgas_cur[igas_h2so4]     = qgas_cur[igas_h2so4] - tmp_q_del;
  }
}  // end mam_newnuc_1subarea
//--------------------------------------------------------------------------------
// Call aerosol microphysics processes for a single (cloudy or clear) subarea
//
// qgas3, qaer3, qaercw3, qnum3, qnumcw3 are the current incoming TMRs
// qgas_cur, qaer_cur, qaercw_cur, qnum_cur, qnumcw_cur are the updated
// outgoing TMRs
//
// In a clear subarea, calculate
//  - gas-aerosol exchange (condensation/evaporation)
//  - growth from smaller to larger modes (renaming) due to condensation
//  - new particle nucleation
//  - coagulation
//  - transfer of particles from hydrophobic modes to hydrophilic modes
//  (aging)
//    due to condensation and coagulation
//
// In a cloudy subarea,
//  - when do_cond = false, this routine only calculate changes involving
//    growth from smaller to larger modes (renaming) following cloud chemistry
//    so gas TMRs are not changed
//  - when do_cond = true, this routine also calculates changes involving
//    gas-aerosol exchange (condensation/evaporation)
//  - transfer of particles from hydrophobic modes to hydrophilic modes
//  (aging)
//       due to condensation
// Currently, in a cloudy subarea, this routine does not do
//  - new particle nucleation - because h2so4 gas conc. should be very low in
//  cloudy air
//  - coagulation - because cloud-borne aerosol would need to be included
//--------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
void mam_amicphys_1subarea(
    // in
    const int newnuc_h2so4_conc_optaa, const int gaexch_h2so4_uptake_optaa,
    const bool do_cond_sub, const bool do_rename_sub, const bool do_newnuc_sub,
    const bool do_coag_sub, const Real deltat, const int jsubarea,
    const bool iscldy_subarea, const Real afracsub, const Real temp,
    const Real pmid, const Real pdel, const Real zmid, const Real pblh,
    const Real relhumsub, const Real (&dgn_a)[num_modes],
    const Real (&dgn_awet)[num_modes], const Real (&wetdens)[num_modes],
    const Real (&qgas1)[max_gas()], const Real (&qgas3)[max_gas()],
    // inout
    Real (&qgas_cur)[max_gas()], Real (&qgas_delaa)[max_gas()][nqtendaa()],
    // in
    const Real (&qnum3)[num_modes],
    // inout
    Real (&qnum_cur)[num_modes], Real (&qnum_delaa)[num_modes][nqtendaa()],
    // in
    const Real (&qaer2)[num_aerosol_ids][num_modes],
    const Real (&qaer3)[num_aerosol_ids][num_modes],
    // inout
    Real (&qaer_cur)[num_aerosol_ids][num_modes],
    Real (&qaer_delaa)[num_aerosol_ids][num_modes][nqtendaa()],
    // in
    const Real (&qwtr3)[num_modes],
    // inout
    Real (&qwtr_cur)[num_modes],
    // in
    const Real (&qnumcw3)[num_modes],
    // inout
    Real (&qnumcw_cur)[num_modes],
    Real (&qnumcw_delaa)[num_modes][nqqcwtendaa()],
    // in
    const Real (&qaercw2)[num_aerosol_ids][num_modes],
    const Real (&qaercw3)[num_aerosol_ids][num_modes],
    // inout
    Real (&qaercw_cur)[num_aerosol_ids][num_modes],
    Real (&qaercw_delaa)[num_aerosol_ids][num_modes][nqqcwtendaa()])

{
  // do_cond_sub, do_rename_sub: true if the aero. microp. process is
  //                             calculated in this subarea
  // do_newnuc_sub, do_coag_sub: true if the aero.  microp. process is
  //                             calculated in this subarea
  // iscldy_subarea: true if sub-area is cloudy
  // kk: level indices
  // jsubarea, nsubarea: sub-area index, number of sub-areas
  // afracsub: fractional area of subarea [unitless, 0-1]
  // deltat: time step [s]
  // temp: air temperature at model levels [K]
  // pmid: air pressure at layer center [Pa]
  // pdel: pressure thickness of layer [Pa]
  // zmid: altitude (above ground) at layer center [m]
  // pblh: planetary boundary layer depth [m]
  // relhum: relative humidity [unitless, 0-1]
  // dgn_a   (max_mode): dry geo. mean diameter [m] of number distribution
  // dgn_awet(max_mode): wet geo. mean diameter [m] of number distribution
  // wetdens (max_mode): interstitial aerosol wet density [kg/m3]

  // Subare mixing ratios qXXXN (X=gas,aer,wat,num; N=1:4):
  //
  //    XXX=gas - gas species [kmol/kmol]
  //    XXX=aer - aerosol mass species (excluding water) [kmol/kmol]
  //    XXX=wat - aerosol water [kmol/kmol]
  //    XXX=num - aerosol number [#/kmol]
  //
  //    N=1 - before gas-phase chemistry
  //    N=2 - before cloud chemistry
  //    N=3 - current incoming values (before gas-aerosol exchange, newnuc,
  //    coag) N=_cur - updated outgoing values (after  gas-aerosol exchange,
  //    newnuc, coag)
  //
  // qgas1, qgas3 [kmol/kmol]
  // qgas_cur     [kmol/kmol]

  // qnum3    [#/kmol]
  // qnum_cur [#/kmol]

  // qaer2, qaer3 [kmol/kmol]
  // qaer_cur[kmol/kmol]

  // qnumcw3[#/kmol]
  // qnumcw_cur  [#/kmol]

  // qaercw2, qaercw3 [kmol/kmol]
  // qaercw_cur  [kmol/kmol]

  // qwtr3      [kmol/kmol]
  // qwtr_cur   [kmol/kmol]

  // qXXX_delaa are TMR changes (increments, not tendencies) of different
  // microphysics processes. These are diagnostics sent to history output;
  // they do not directly affect time integration.

  // qgas_delaa   [kmol/kmol]
  // qnum_delaa   [   #/kmol]
  // qaer_delaa   [kmol/kmol]
  // qnumcw_delaa [   #/kmol]
  // qaercw_delaa [kmol/kmol]

  // type ( misc_vars_aa_type ), intent(inout) :: misc_vars_aa_sub

  //---------------------------------------------------------------------------------------
  // Calculate air molar density [kmol/m3] to be passed on to individual
  // parameterizations
  //---------------------------------------------------------------------------------------
  // BAD CONSTANT
  // Universal gas constant (J/K/kmol)
  constexpr Real r_universal = 8314.46759100000;
  const Real aircon          = pmid / (r_universal * temp);

  //----------------------------------------------------------
  // Initializ mixing ratios with the before-amicphys values
  //----------------------------------------------------------

  copy_1d_array(max_gas(), qgas3,  // in
                qgas_cur);         // out

  constexpr int nspecies = num_aerosol_ids;
  constexpr int nmodes   = num_modes;

  copy_2d_array(nspecies, nmodes, qaer3,  // in
                qaer_cur);                // out

  copy_1d_array(nmodes, qnum3,  // in
                qnum_cur);      // out

  copy_1d_array(nmodes, qwtr3,  // in
                qwtr_cur);      // out

  if(iscldy_subarea) {
    copy_1d_array(nmodes, qnumcw3,            // in
                  qnumcw_cur);                // out
    copy_2d_array(nspecies, nmodes, qaercw3,  // in
                  qaercw_cur);                // out
  }                                           // iscldy_subarea

  //---------------------------------------------------------------------
  // Diagnose net production rate of H2SO4 gas production
  // cause by other processes (e.g., gas chemistry and cloud chemistry)
  //---------------------------------------------------------------------
  Real qgas_netprod_otrproc[max_gas()] = {0};
  assign_1d_array(max_gas(), 0.0,         // in
                  qgas_netprod_otrproc);  // out

  // If gaexch_h2so4_uptake_optaa == 2, then
  //  - if qgas increases from pre-gaschem to post-cldchem,
  //    start from the pre-gaschem mix-ratio and add in the production during
  //    the integration
  //  - if it decreases,  start from post-cldchem mix-ratio

  if((do_cond_sub) && (gaexch_h2so4_uptake_optaa == 2)) {
    for(int igas = 0; igas < max_gas(); ++igas) {
      if((igas == igas_h2so4) || (igas == igas_nh3)) {
        qgas_netprod_otrproc[igas] = (qgas3[igas] - qgas1[igas]) / deltat;
        qgas_cur[igas] = (qgas_netprod_otrproc[igas] >= 0) ? qgas1[igas] : 0;
      }  // h2so4, igas_nh3
    }    // igas
  }      // do_cond_sub,gaexch_h2so4_uptake_optaa

  constexpr int ntsubstep = 1;
  const Real del_h2so4_gasprod =
      haero::max(qgas3[igas_h2so4] - qgas1[igas_h2so4], 0) / ntsubstep;
  //-----------------------------------
  // Initialize increment diagnostics
  //-----------------------------------

  assign_2d_array(max_gas(), nqtendaa(), 0,  // in
                  qgas_delaa);               // out

  assign_2d_array(nmodes, nqtendaa(), 0,  // in
                  qnum_delaa);            // out

  assign_3d_array(nspecies, nmodes, nqtendaa(), 0,  // in
                  qaer_delaa);                      // out

  assign_2d_array(nmodes, nqqcwtendaa(), 0,  // in
                  qnumcw_delaa);             // out

  assign_3d_array(nspecies, nmodes, nqqcwtendaa(), 0,  // in
                  qaercw_delaa);                       // out

  Real ncluster_tend_nnuc_1grid = 0;

  //***********************************
  // loop over multiple time sub-steps
  //***********************************
  EKAT_KERNEL_ASSERT_MSG(ntsubstep != 0,
                         "Error! mam_amicphys_1subarea: "
                         " ntsubstep should not be equal to 0\n");
  const int dtsubstep = deltat / ntsubstep;

  Real qgas_sv1[max_gas()];
  Real qnum_sv1[nmodes];
  Real qaer_sv1[nspecies][nmodes];

  Real del_h2so4_aeruptk;    // [kmol/kmol]
  Real qgas_avg[max_gas()];  // [kmol/kmol]

  // Mixing ratio increments of sub-timesteps used for process coupling

  Real qnum_delsub_cond[nmodes];                   // [   #/kmol]
  Real qnum_delsub_coag[nmodes];                   // [   #/kmol]
  Real qaer_delsub_cond[nspecies][nmodes];         // [   #/kmol]
  Real qaer_delsub_coag[nspecies][nmodes];         // [kmol/kmol]
  Real qaer_delsub_grow4rnam[nspecies][nmodes];    // [kmol/kmol]
  Real qaercw_delsub_grow4rnam[nspecies][nmodes];  // [kmol/kmol]

  constexpr int max_agepair = AeroConfig::max_agepair();
  Real qaer_delsub_coag_in[nspecies][max_agepair];  // [kmol/kmol]

  for(int jtsubstep = 0; jtsubstep < ntsubstep; ++jtsubstep) {
    //======================
    // Gas-aerosol exchange
    //======================
    Real uptkrate_h2so4 = 0;

    if(do_cond_sub) {
      copy_1d_array(max_gas(), qgas_cur,  // in
                    qgas_sv1);            // out
      copy_1d_array(nmodes, qnum_cur,     // in
                    qnum_sv1);            // out

      copy_2d_array(nspecies, nmodes, qaer_cur,  // in
                    qaer_sv1);                   // out

      // max_mode in MAM4 is different from nmodes(max_mode = nmodes+1)
      // Here we create temporary arrays for now, but we should make
      // it consistent to avoid this extra memoery
      constexpr int max_mode = nmodes + 1;
      Real qaer_cur_tmp[nspecies][max_mode];
      Real qnum_cur_tmp[max_mode];
      Real qwtr_cur_tmp[max_mode];

      // NOTE: we cannot use copy_2d_array here as arrays extent is max_mode
      //  but we are copying till nmodes
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_cur_tmp[is][im] = qaer_cur[is][im];
        }
      }

      copy_1d_array(nmodes, qwtr_cur,  // in
                    qwtr_cur_tmp);     // out

      copy_1d_array(nmodes, qnum_cur,  // in
                    qnum_cur_tmp);     // out

      Real uptkaer[max_gas()][max_mode];

      mam4::mam_gasaerexch_1subarea(
          jtsubstep, dtsubstep, temp, pmid, aircon, nmodes,  // in
          qgas_cur, qgas_avg,                                // inout
          qgas_netprod_otrproc,                              // in
          qaer_cur_tmp, qnum_cur_tmp, qwtr_cur_tmp,          // inout
          dgn_awet,                                          // in
          uptkaer, uptkrate_h2so4);                          // inout

      // copy back the values for aerosols
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_cur[is][im] = qaer_cur_tmp[is][im];
        }
      }

      copy_1d_array(nmodes, qwtr_cur_tmp,  // in
                    qwtr_cur);             // out

      copy_1d_array(nmodes, qnum_cur_tmp,  // in
                    qnum_cur);             // out

      for(int ig = 0; ig < max_gas(); ++ig) {
        qgas_delaa[ig][iqtend_cond()] =
            qgas_delaa[ig][iqtend_cond()] +
            (qgas_cur[ig] -
             (qgas_sv1[ig] + qgas_netprod_otrproc[ig] * dtsubstep));
      }
      for(int im = 0; im < nmodes; ++im) {
        qnum_delsub_cond[im] = qnum_cur[im] - qnum_sv1[im];
      }
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_delsub_cond[is][im] = qaer_cur[is][im] - qaer_sv1[is][im];
        }
      }

      del_h2so4_aeruptk =
          qgas_cur[igas_h2so4] -
          (qgas_sv1[igas_h2so4] + qgas_netprod_otrproc[igas_h2so4] * dtsubstep);

    } else {                              // do_cond_sub
      copy_1d_array(max_gas(), qgas_cur,  // in
                    qgas_avg);            // out

      assign_2d_array(nspecies, nmodes, 0,  // in
                      qaer_delsub_cond);    // out

      assign_1d_array(nmodes, 0.0,        // in
                      qnum_delsub_cond);  // out
      del_h2so4_aeruptk = 0;

    }  // do_cond_sub

    //====================================
    // Renaming after "continuous growth"
    //====================================
    if(do_rename_sub) {
      constexpr int dest_mode_of_mode[nmodes] = {-1, 0, -1, -1};

      //---------------------------------------------------------
      // Calculate changes in aerosol mass mixing ratios due to
      //  - gas condensation/evaporation
      //  - cloud chemistry (if the subarea is cloudy)
      //---------------------------------------------------------
      copy_2d_array(nspecies, nmodes, qaer_delsub_cond,  // in
                    qaer_delsub_grow4rnam);              // out

      if(iscldy_subarea) {
        for(int is = 0; is < nspecies; ++is) {
          for(int im = 0; im < nmodes; ++im) {
            qaer_delsub_grow4rnam[is][im] =
                (qaer3[is][im] - qaer2[is][im]) / ntsubstep +
                qaer_delsub_grow4rnam[is][im];
            qaercw_delsub_grow4rnam[is][im] =
                (qaercw3[is][im] - qaercw2[is][im]) / ntsubstep;
          }
        }
      }

      //----------
      // Renaming
      //----------
      copy_1d_array(nmodes, qnum_cur,  // in
                    qnum_sv1);         // out

      copy_2d_array(nspecies, nmodes, qaer_cur,  // in
                    qaer_sv1);                   // out

      Real qnumcw_sv1[nmodes];
      copy_1d_array(nmodes, qnumcw_cur,  // in
                    qnumcw_sv1);         // out
      Real qaercw_sv1[nspecies][nmodes];
      copy_2d_array(nspecies, nmodes, qaercw_cur,  // in
                    qaercw_sv1);                   // out

      Real mean_std_dev[nmodes];
      Real fmode_dist_tail_fac[nmodes];
      Real v2n_lo_rlx[nmodes];
      Real v2n_hi_rlx[nmodes];
      Real ln_diameter_tail_fac[nmodes];
      int num_pairs = 0;
      Real diameter_cutoff[nmodes];
      Real ln_dia_cutoff[nmodes];
      Real diameter_threshold[nmodes];

      rename::find_renaming_pairs(
          dest_mode_of_mode,                              // in
          mean_std_dev, fmode_dist_tail_fac, v2n_lo_rlx,  // out
          v2n_hi_rlx, ln_diameter_tail_fac, num_pairs,    // out
          diameter_cutoff, ln_dia_cutoff,                 // out
          diameter_threshold);                            // out
      Real dgnum_amode[nmodes];
      for(int m = 0; m < nmodes; ++m) {
        dgnum_amode[m] = modes(m).nom_diameter;
      }
      // BAD_CONSTANT
      constexpr Real smallest_dryvol_value = 1.0e-25;

      // swap dimensions as mam_rename_1subarea_ uses output arrays in
      //  a swapped dimension order
      Real qaer_cur_tmp[nmodes][nspecies];
      Real qaer_delsub_grow4rnam_tmp[nmodes][nspecies];
      Real qaercw_cur_tmp[nmodes][nspecies];
      Real qaercw_delsub_grow4rnam_tmp[nmodes][nspecies];
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_cur_tmp[im][is]                = qaer_cur[is][im];
          qaer_delsub_grow4rnam_tmp[im][is]   = qaer_delsub_grow4rnam[is][im];
          qaercw_cur_tmp[im][is]              = qaercw_cur[is][im];
          qaercw_delsub_grow4rnam_tmp[im][is] = qaercw_delsub_grow4rnam[is][im];
        }
      }
      Rename rename;
      rename.mam_rename_1subarea_(
          iscldy_subarea, smallest_dryvol_value, dest_mode_of_mode,    // in
          mean_std_dev, fmode_dist_tail_fac, v2n_lo_rlx, v2n_hi_rlx,   // in
          ln_diameter_tail_fac, num_pairs, diameter_cutoff,            // in
          ln_dia_cutoff, diameter_threshold, mass_2_vol, dgnum_amode,  // in
          qnum_cur, qaer_cur_tmp,                                      // out
          qaer_delsub_grow4rnam_tmp,                                   // in
          qnumcw_cur, qaercw_cur_tmp,                                  // out
          qaercw_delsub_grow4rnam_tmp);                                // in

      // copy the output back to the variables
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_cur[is][im]                = qaer_cur_tmp[im][is];
          qaer_delsub_grow4rnam[is][im]   = qaer_delsub_grow4rnam_tmp[im][is];
          qaercw_cur[is][im]              = qaercw_cur_tmp[im][is];
          qaercw_delsub_grow4rnam[is][im] = qaercw_delsub_grow4rnam_tmp[im][is];
        }
      }

      //------------------------
      // Accumulate increments
      //------------------------
      for(int im = 0; im < nmodes; ++im) {
        qnum_delaa[im][iqtend_rnam()] =
            qnum_delaa[im][iqtend_rnam()] + (qnum_cur[im] - qnum_sv1[im]);
      }

      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_delaa[is][im][iqtend_rnam()] =
              qaer_delaa[is][im][iqtend_rnam()] +
              (qaer_cur[is][im] - qaer_sv1[is][im]);
        }
      }

      if(iscldy_subarea) {
        for(int im = 0; im < nmodes; ++im) {
          qnumcw_delaa[im][iqqcwtend_rnam()] =
              qnumcw_delaa[im][iqqcwtend_rnam()] +
              (qnumcw_cur[im] - qnumcw_sv1[im]);
        }
      }  // if iscldy_subarea

      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaercw_delaa[is][im][iqqcwtend_rnam()] =
              qaercw_delaa[is][im][iqqcwtend_rnam()] +
              (qaercw_cur[is][im] - qaercw_sv1[is][im]);
        }
      }
    }  // do_rename_sub

    //====================================
    // New particle formation (nucleation)
    //====================================
    if(do_newnuc_sub) {
      copy_1d_array(max_gas(), qgas_cur,  // in
                    qgas_sv1);            // out
      copy_1d_array(nmodes, qnum_cur,     // in
                    qnum_sv1);            // out

      copy_2d_array(nspecies, nmodes, qaer_cur,  // in
                    qaer_sv1);                   // out

      Real dnclusterdt_substep;
      mam_newnuc_1subarea(igas_h2so4, gaexch_h2so4_uptake_optaa,
                          newnuc_h2so4_conc_optaa, jsubarea, dtsubstep,  // in
                          temp,                                          // in
                          pmid, aircon, zmid, pblh,                      // in
                          relhumsub, uptkrate_h2so4, del_h2so4_gasprod,  // in
                          del_h2so4_aeruptk,                             // in
                          qgas_cur, qgas_avg, qnum_cur, qaer_cur,        // out
                          dnclusterdt_substep);                          // out

      for(int ig = 0; ig < max_gas(); ++ig) {
        qgas_delaa[ig][iqtend_nnuc()] += (qgas_cur[ig] - qgas_sv1[ig]);
      }
      for(int im = 0; im < nmodes; ++im) {
        qnum_delaa[im][iqtend_nnuc()] += (qnum_cur[im] - qnum_sv1[im]);
      }
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_delaa[is][im][iqtend_nnuc()] +=
              (qaer_cur[is][im] - qaer_sv1[is][im]);
        }
      }
      EKAT_KERNEL_ASSERT_MSG(deltat != 0,
                             "Error! mam_amicphys_1subarea: "
                             "deltat should not be equal to zero \n");
      ncluster_tend_nnuc_1grid += dnclusterdt_substep * (dtsubstep / deltat);

    }  // do_newnuc_sub

    //====================================
    // Coagulation
    //====================================
    if(do_coag_sub) {
      copy_1d_array(nmodes, qnum_cur,  // in
                    qnum_sv1);         // out

      copy_2d_array(nspecies, nmodes, qaer_cur,  // in
                    qaer_sv1);                   // out

      mam4::coagulation::mam_coag_1subarea(
          dtsubstep, temp, pmid, aircon,             // in
          dgn_awet, wetdens,                         // in
          qnum_cur, qaer_cur, qaer_delsub_coag_in);  // inout, inout, out

      for(int im = 0; im < nmodes; ++im) {
        qnum_delsub_coag[im] = qnum_cur[im] - qnum_sv1[im];
      }

      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_delsub_coag[is][im] = qaer_cur[is][im] - qaer_sv1[is][im];
        }
      }

      for(int im = 0; im < nmodes; ++im) {
        qnum_delaa[im][iqtend_coag()] += qnum_delsub_coag[im];
      }
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_delaa[is][im][iqtend_coag()] += qaer_delsub_coag[is][im];
        }
      }

    } else {
      assign_2d_array(nspecies, max_agepair, 0.0,  // in
                      qaer_delsub_coag_in);        // out
      assign_2d_array(nspecies, nmodes, 0.0,       // in
                      qaer_delsub_coag);           // out
      assign_1d_array(nmodes, 0.0,                 // in
                      qnum_delsub_coag);           // out

    }  // do_coag_sub

    //====================================
    // primary carbon aging
    //====================================
    const bool do_aging_in_subarea =
        (n_agepair() > 0) &&
        ((!iscldy_subarea) || (iscldy_subarea && do_cond_sub));

    if(do_aging_in_subarea) {
      mam4::aging::mam_pcarbon_aging_1subarea(
          dgn_a,                                         // input
          qnum_cur, qnum_delsub_cond, qnum_delsub_coag,  // in-outs
          qaer_cur, qaer_delsub_cond, qaer_delsub_coag,  // in-outs
          qaer_delsub_coag_in);                          // in-outs
    }                                                    // do_aging_in_subarea

    // The following block has to be placed here (after both condensation and
    // aging) as both can change the values of qnum_delsub_cond and
    // qaer_delsub_cond.

    if(do_cond_sub) {
      for(int im = 0; im < nmodes; ++im) {
        qnum_delaa[im][iqtend_cond()] =
            qnum_delaa[im][iqtend_cond()] + qnum_delsub_cond[im];
      }
      for(int is = 0; is < nspecies; ++is) {
        for(int im = 0; im < nmodes; ++im) {
          qaer_delaa[is][im][iqtend_cond()] =
              qaer_delaa[is][im][iqtend_cond()] + qaer_delsub_cond[is][im];
        }
      }
    }  // do_cond_sub

  }  // jtsubstep_loop

}  // mam_amicphys_1subarea

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void mam_amicphys_1gridcell(
    // in
    const AmicPhysConfig &config, const Real deltat, const int nsubarea,
    const int ncldy_subarea, const bool (&iscldy_subarea)[maxsubarea()],
    const Real (&afracsub)[maxsubarea()], const Real temp, const Real pmid,
    const Real pdel, const Real zmid, const Real pblh,
    const Real (&relhumsub)[maxsubarea()], const Real (&dgn_a)[num_modes],
    const Real (&dgn_awet)[num_modes], const Real (&wetdens)[num_modes],
    const Real (&qsub1)[gas_pcnst][maxsubarea()],
    const Real (&qsub2)[gas_pcnst][maxsubarea()],
    const Real (&qqcwsub2)[gas_pcnst][maxsubarea()],
    const Real (&qsub3)[gas_pcnst][maxsubarea()],
    const Real (&qqcwsub3)[gas_pcnst][maxsubarea()],
    const Real (&qaerwatsub3)[num_modes][maxsubarea()],
    // out
    Real (&qsub4)[gas_pcnst][maxsubarea()],
    Real (&qqcwsub4)[gas_pcnst][maxsubarea()],
    Real (&qaerwatsub4)[num_modes][maxsubarea()],
    Real (&qsub_tendaa)[gas_pcnst][nqtendaa()][maxsubarea()],
    Real (&qqcwsub_tendaa)[gas_pcnst][nqqcwtendaa()][maxsubarea()]) {
  //
  // calculates changes to gas and aerosol sub-area TMRs (tracer mixing
  // ratios) qsub3 and qqcwsub3 are the incoming current TMRs qsub4 and
  // qqcwsub4 are the outgoing updated TMRs
  //
  // qsubN and qqcwsubN (N=1:4) are tracer mixing ratios (TMRs, mol/mol or
  // #/kmol) in sub-areas
  //    currently there are just clear and cloudy sub-areas
  //    the N=1:4 have same meanings as for qgcmN
  //    N=1 - before gas-phase chemistry
  //    N=2 - before cloud chemistry
  //    N=3 - incoming values (before gas-aerosol exchange, newnuc, coag)
  //    N=4 - outgoing values (after  gas-aerosol exchange, newnuc, coag)
  // qsub_tendaa and qqcwsub_tendaa are TMR tendencies
  //    for different processes, which are used to produce history output
  // the processes are condensation/evaporation (and associated aging),
  //    renaming, coagulation, and nucleation

  constexpr int mdo_gaexch_cldy_subarea = 0;
  // the qq--4 values will be equal to qq--3 values unless they get changed
  for(int i = 0; i < num_gas_ids; ++i) {
    for(int j = 1; j < maxsubarea(); ++j) {
      qsub4[i][j]    = qsub3[i][j];
      qqcwsub4[i][j] = qqcwsub3[i][j];
    }
  }

  for(int i = 0; i < num_modes; ++i) {
    for(int j = 0; j < maxsubarea(); ++j) {
      qaerwatsub4[i][j] = qaerwatsub3[i][j];
    }
  }

  assign_3d_array(num_gas_ids, nqtendaa(), maxsubarea(), 0.0,  // in
                  qsub_tendaa);                                // out

  assign_3d_array(num_gas_ids, nqqcwtendaa(), maxsubarea(), 0.0,  // in
                  qqcwsub_tendaa);                                // out

  EKAT_KERNEL_ASSERT_MSG(nsubarea < maxsubarea(),
                         "Error! mam_amicphys_1gridcell: "
                         "nsubarea should be < maxsubarea() \n");
  for(int jsub = 1; jsub <= nsubarea; ++jsub) {
    bool do_cond;
    bool do_rename;
    bool do_newnuc;
    bool do_coag;
    if(iscldy_subarea[jsub]) {
      do_cond   = config.do_cond;
      do_rename = config.do_rename;
      do_newnuc = false;
      do_coag   = false;
      if(mdo_gaexch_cldy_subarea <= 0) do_cond = false;
    } else {
      do_cond   = config.do_cond;
      do_rename = config.do_rename;
      do_newnuc = config.do_newnuc;
      do_coag   = config.do_coag;
    }
    const bool do_map_gas_sub = do_cond || do_newnuc;

    // map incoming sub-area mix-ratios to gas/aer/num arrays
    Real qgas1[max_gas()] = {0};
    Real qgas2[max_gas()] = {0};
    Real qgas3[max_gas()] = {0};
    Real qgas4[max_gas()] = {0};
    assign_1d_array(max_gas(), 0.0,  // in
                    qgas1);          // out
    assign_1d_array(max_gas(), 0.0,  // in
                    qgas2);          // out
    assign_1d_array(max_gas(), 0.0,  // in
                    qgas3);          // out
    assign_1d_array(max_gas(), 0.0,  // in
                    qgas4);          // out

    if(do_map_gas_sub) {
      // for cldy subarea, only do gases if doing gaexch
      for(int igas = 0; igas < max_gas(); ++igas) {
        const int l = lmap_gas(igas);
        qgas1[igas] = qsub1[l][jsub] * fcvt_gas(igas);
        qgas2[igas] = qsub2[l][jsub] * fcvt_gas(igas);
        qgas3[igas] = qsub3[l][jsub] * fcvt_gas(igas);
        qgas4[igas] = qgas3[igas];
      }
    }

    Real qaer2[num_aerosol_ids][num_modes] = {0};
    Real qnum2[num_modes]                  = {0};
    Real qaer3[num_aerosol_ids][num_modes] = {0};
    Real qnum3[num_modes]                  = {0};
    Real qaer4[num_aerosol_ids][num_modes] = {0};
    Real qnum4[num_modes]                  = {0};
    Real qwtr3[num_modes]                  = {0};
    Real qwtr4[num_modes]                  = {0};

    assign_2d_array(num_aerosol_ids, num_modes, 0.0,  // in
                    qaer2);                           // out
    assign_2d_array(num_aerosol_ids, num_modes, 0.0,  // in
                    qaer3);                           // out
    assign_2d_array(num_aerosol_ids, num_modes, 0.0,  // in
                    qaer4);                           // out

    assign_1d_array(num_modes, 0.0,  // in
                    qnum2);          // out
    assign_1d_array(num_modes, 0.0,  // in
                    qnum3);          // out
    assign_1d_array(num_modes, 0.0,  // in
                    qnum4);          // out
    assign_1d_array(num_modes, 0.0,  // in
                    qwtr3);          // out
    assign_1d_array(num_modes, 0.0,  // in
                    qwtr4);          // out

    for(int imode = 0; imode < num_modes; ++imode) {
      const int ln = lmap_num(imode);
      qnum2[imode] = qsub2[ln][jsub] * fcvt_num();
      qnum3[imode] = qsub3[ln][jsub] * fcvt_num();
      qnum4[imode] = qnum3[imode];
      for(int iaer = 0; iaer < num_aerosol_ids; ++iaer) {
        const int la = lmap_aer(iaer, imode);
        if(la > 0) {
          qaer2[iaer][imode] = qsub2[la][jsub] * fcvt_aer(iaer);
          qaer3[iaer][imode] = qsub3[la][jsub] * fcvt_aer(iaer);
          qaer4[iaer][imode] = qaer3[iaer][imode];
        }
      }  // for iaer
      qwtr3[imode] = qaerwatsub3[imode][jsub] * fcvt_wtr();
      qwtr4[imode] = qwtr3[imode];
    }  // for imode

    Real qaercw2[num_aerosol_ids][num_modes] = {0};
    Real qnumcw2[num_modes]                  = {0};
    Real qaercw3[num_aerosol_ids][num_modes] = {0};
    Real qnumcw3[num_modes]                  = {0};
    Real qaercw4[num_aerosol_ids][num_modes] = {0};
    Real qnumcw4[num_modes]                  = {0};

    assign_2d_array(num_aerosol_ids, num_modes, 0.0,  // in
                    qaercw2);                         // out
    assign_2d_array(num_aerosol_ids, num_modes, 0.0,  // in
                    qaercw3);                         // out
    assign_2d_array(num_aerosol_ids, num_modes, 0.0,  // in
                    qaercw4);                         // out

    assign_1d_array(num_modes, 0.0,  // in
                    qnumcw2);        // out
    assign_1d_array(num_modes, 0.0,  // in
                    qnumcw3);        // out
    assign_1d_array(num_modes, 0.0,  // in
                    qnumcw4);        // out

    if(iscldy_subarea[jsub]) {
      for(int imode = 0; imode < num_modes; ++imode) {
        qnumcw2[imode] = 0;
        qnumcw3[imode] = 0;
        qnumcw4[imode] = 0;
      }  // imode
      for(int iaer = 0; iaer < num_aerosol_ids; ++iaer) {
        for(int imode = 0; imode < num_modes; ++imode) {
          qaercw2[iaer][imode] = 0;
          qaercw3[iaer][imode] = 0;
          qaercw4[iaer][imode] = 0;
        }  // imode
      }    // iaer
      // only do cloud-borne for cloudy
      for(int imode = 0; imode < num_modes; ++imode) {
        int ln         = lmap_numcw(imode);
        qnumcw2[imode] = qqcwsub2[ln][jsub] * fcvt_num();
        qnumcw3[imode] = qqcwsub3[ln][jsub] * fcvt_num();
        qnumcw4[imode] = qnumcw3[imode];
      }  // imode
      for(int iaer = 0; iaer < num_aerosol_ids; ++iaer) {
        for(int imode = 0; imode < num_modes; ++imode) {
          int la = lmap_aer(iaer, imode);
          if(la > 0) {
            qaercw2[iaer][imode] = qqcwsub2[la][jsub] * fcvt_aer(iaer);
            qaercw3[iaer][imode] = qqcwsub3[la][jsub] * fcvt_aer(iaer);
            qaercw4[iaer][imode] = qaercw3[iaer][imode];
          }  // la
        }    // imode
      }      // iaer
    }        // iscldy_subarea

    Real qgas_delaa[max_gas()][nqtendaa()]                       = {};
    Real qnum_delaa[num_modes][nqtendaa()]                       = {};
    Real qnumcw_delaa[num_modes][nqqcwtendaa()]                  = {};
    Real qaer_delaa[num_aerosol_ids][num_modes][nqtendaa()]      = {};
    Real qaercw_delaa[num_aerosol_ids][num_modes][nqqcwtendaa()] = {};

    mam_amicphys_1subarea(
        // in
        config.gaexch_h2so4_uptake_optaa, config.newnuc_h2so4_conc_optaa,
        do_cond, do_rename, do_newnuc, do_coag, deltat, jsub,
        iscldy_subarea[jsub], afracsub[jsub], temp, pmid, pdel, zmid, pblh,
        relhumsub[jsub], dgn_a, dgn_awet, wetdens, qgas1, qgas3, qgas4,
        qgas_delaa,              // out
        qnum3,                   // in
        qnum4, qnum_delaa,       // out
        qaer2, qaer3,            // in
        qaer4, qaer_delaa,       // out
        qwtr3,                   // in
        qwtr4,                   // out
        qnumcw3,                 // in
        qnumcw4, qnumcw_delaa,   // out
        qaercw2, qaercw3,        // in
        qaercw4, qaercw_delaa);  // out

    // FIXME: Enable this functionality
    /*if (nsubarea == 1 || !iscldy_subarea[jsub]) {
     ncluster_tend_nnuc_1grid = ncluster_tend_nnuc_1grid &
                                           +
    misc_vars_aa_sub(jsub)%ncluster_tend_nnuc_1grid*afracsub(jsub)
    }*/

    // map gas/aer/num arrays (mix-ratio and del=change) back to sub-area arrays

    if(do_map_gas_sub) {
      for(int igas = 0; igas < max_gas(); ++igas) {
        int ll          = lmap_gas(igas);
        qsub4[ll][jsub] = qgas4[igas] / fcvt_gas(igas);
        for(int jj = 0; jj < nqtendaa(); ++jj) {
          qsub_tendaa[ll][jj][jsub] =
              qgas_delaa[igas][jj] / (fcvt_gas(igas) * deltat);
        }
      }  // igas
    }    // do_map_gas_sub

    for(int imode = 0; imode < num_modes; ++imode) {
      int ll          = lmap_num(imode);
      qsub4[ll][jsub] = qnum4[imode] / (fcvt_num());
      for(int jj = 0; jj < nqtendaa(); ++jj) {
        qsub_tendaa[ll][jj][jsub] =
            qnum_delaa[imode][jj] / (fcvt_num() * deltat);
      }
      for(int iaer = 0; iaer < num_aerosol_ids; ++iaer) {
        int la = lmap_aer(iaer, imode);
        if(la > 0) {
          qsub4[la][jsub] = qaer4[iaer][imode] / fcvt_aer(iaer);
          for(int jj = 0; jj < nqtendaa(); ++jj) {
            qsub_tendaa[la][jj][jsub] =
                qaer_delaa[iaer][imode][jj] / (fcvt_aer(iaer) * deltat);
          }  // jj
        }    // la
      }      // iaer
      qaerwatsub4[imode][jsub] = qwtr4[imode] / fcvt_wtr();

      if(iscldy_subarea[jsub]) {
        int lc             = lmap_numcw(imode);
        qqcwsub4[lc][jsub] = qnumcw4[imode] / fcvt_num();
        for(int jj = 0; jj < nqqcwtendaa(); ++jj) {
          qqcwsub_tendaa[lc][jj][jsub] =
              qnumcw_delaa[imode][jj] / (fcvt_num() * deltat);
        }  // jj
        for(int iaer = 0; iaer < num_aerosol_ids; ++iaer) {
          int lca = lmap_aercw(iaer, imode);
          if(lca > 0) {
            qqcwsub4[lca][jsub] = qaercw4[iaer][imode] / fcvt_aer(iaer);
            for(int jj = 0; jj < nqqcwtendaa(); ++jj) {
              qqcwsub_tendaa[lca][jj][jsub] =
                  qaercw_delaa[iaer][imode][jj] / (fcvt_aer(iaer) * deltat);
            }  // jj
          }    // lca
        }      // iaer
      }        // iscldy_subarea
    }          // imode
  }            // main_jsub_loop

}  // mam_amicphys_1gridcell

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void form_gcm_of_gases_and_aerosols_from_subareas(
    // in
    const int nsubarea, const int ncldy_subarea,
    const Real (&afracsub)[maxsubarea()],
    const Real (&qsub)[gas_pcnst][maxsubarea()],
    const Real (&qqcwsub)[gas_pcnst][maxsubarea()],
    const Real (&qqcwgcm_old)[gas_pcnst],
    // out
    Real (&qgcm)[gas_pcnst], Real (&qqcwgcm)[gas_pcnst]) {
  //--------------------------------------------------------------------------
  // Purpose: Form grid cell mean values by calculating area-weighted averages
  // of the subareas.
  //           - For gases and interstitial aerosols, sum over all active
  //           subareas.
  //           - For cloud-borne aerosols,
  //---------------------------------------------------------------------------

  // nsubarea: # of active subareas
  // ncldy_subarea : # of cloudy subareas
  // afracsub(maxsubarea):area fraction of subareas [unitless]

  // The following arguments are mixing ratios. Their units do not matter for
  // this subroutine.

  // qsub   (ncnst, maxsubarea):gas and interst. aerosol mixing ratios in
  // subareas
  // qqcwsub(ncnst, maxsubarea): cloud-borne aerosol mixing ratios in
  // subareas
  // qqcwgcm_old(ncnst): grid cell mean cloud-borne aerosol mixing
  // ratios before aerosol microphysics calculations
  // qgcm   (ncnst): grid cell mean gas and interst.  aerosol mixing ratios
  // qqcwgcm(ncnst): cloud-borne aerosol mixing ratios in subareas

  // Gases and interstitial aerosols
  assign_1d_array(gas_pcnst, 0.0,  // in
                  qgcm);           // out

  EKAT_KERNEL_ASSERT_MSG(nsubarea < maxsubarea(),
                         "Error! form_gcm_of_gases_and_aerosols_from_subareas: "
                         "nsubarea should be < maxsubarea() \n");

  for(int jsub = 1; jsub <= nsubarea; ++jsub) {
    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      qgcm[icnst] += qsub[icnst][jsub] * afracsub[jsub];
    }
  }

  for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
    qgcm[icnst] = haero::max(0, qgcm[icnst]);
  }

  // Cloud-borne aerosols
  if(ncldy_subarea <= 0) {
    for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
      qqcwgcm[icnst] = qqcwgcm_old[icnst];
    }
  } else {
    assign_1d_array(gas_pcnst, 0.0,  // in
                    qqcwgcm);        // out
    for(int jsub = 1; jsub <= nsubarea; ++jsub) {
      for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
        qqcwgcm[icnst] += qqcwsub[icnst][jsub] * afracsub[jsub];
      }
    }
  }  // if ncldy_subarea
}  // form_gcm_of_gases_and_aerosols_from_subareas

//--------------------------------------------------------------------------------
// Purpose: Get grid cell mean tendencies by calculating area-weighted averages
//          of the values in different subareas.
//--------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
void get_gcm_tend_diags_from_subareas(
    // in
    const int nsubarea, const int ncldy_subarea,
    const Real (&afracsub)[maxsubarea()],
    const Real (&qsub_tendaa)[gas_pcnst][nqtendaa()][maxsubarea()],
    const Real (&qqcwsub_tendaa)[gas_pcnst][nqqcwtendaa()][maxsubarea()],
    // out
    Real (&qgcm_tendaa)[gas_pcnst][nqtendaa()],
    Real (&qqcwgcm_tendaa)[gas_pcnst][nqqcwtendaa()]) {
  // nsubarea: # of active subareas
  // ncldy_subarea: # of cloudy subareas
  // afracsub(maxsubarea): area fraction of subareas [unitless]

  // Gases and interstitial aerosols
  assign_2d_array(gas_pcnst, nqtendaa(), 0.0,  // in
                  qgcm_tendaa);                // out

  EKAT_KERNEL_ASSERT_MSG(nsubarea < maxsubarea(),
                         "Error! get_gcm_tend_diags_from_subareas: "
                         "nsubarea should be < maxsubarea() \n");

  for(int jsub = 1; jsub <= nsubarea; ++jsub) {
    for(int iq = 0; iq < nqtendaa(); ++iq) {
      for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
        qgcm_tendaa[icnst][iq] += qsub_tendaa[icnst][iq][jsub] * afracsub[jsub];
      }
    }
  }

  // Cloud-borne aerosols

  assign_2d_array(gas_pcnst, nqqcwtendaa(), 0.0,  // in
                  qqcwgcm_tendaa);                // out
  if(ncldy_subarea > 0) {
    for(int jsub = 1; jsub <= nsubarea; ++jsub) {
      for(int iq = 0; iq < nqqcwtendaa(); ++iq) {
        for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
          qqcwgcm_tendaa[icnst][iq] +=
              qqcwsub_tendaa[icnst][iq][jsub] * afracsub[jsub];
        }
      }
    }
  }  // if (ncldy_subarea
}  // get_gcm_tend_diags_from_subareas

}  // anonymous namespace

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

KOKKOS_INLINE_FUNCTION
void modal_aero_amicphys_intr(
    // in
    const AmicPhysConfig &config, const Real deltat, const Real temp,
    const Real pmid, const Real pdel, const Real zm, const Real pblh,
    const Real qv, const Real cld,
    // out
    Real qq[gas_pcnst], Real qqcw[gas_pcnst],
    // in
    const Real (&q_pregaschem)[gas_pcnst],
    const Real (&q_precldchem)[gas_pcnst],
    const Real (&qqcw_precldchem)[gas_pcnst], const Real (&dgncur_a)[num_modes],
    const Real (&dgncur_awet)[num_modes],
    const Real (&wetdens_host)[num_modes]) {
  // deltat: time step
  // qq(ncol,pver,pcnst): current tracer mixing ratios (TMRs)
  //                           these values are updated (so out /= in)
  //                        *** MUST BE  #/kmol-air for number
  //                        *** MUST BE mol/mol-air for mass
  //                        *** NOTE ncol dimension
  // qqcw(ncol,pver,pcnst):
  //                          like qq but for cloud-borner tracers
  //                         these values are updated
  // q_pregaschem(ncol,pver,pcnst): qq TMRs    before gas-phase
  // chemistry
  //
  // q_precldchem(ncol,pver,pcnstxx): qq TMRs    before cloud
  // chemistry
  // qqcw_precldchem(ncol,pver,pcnstxx): qqcw TMRs before cloud
  // chemistry
  // t(pcols,pver): temperature at model levels (K)
  // pmid(pcols,pver):pressure at model level centers (Pa)
  // pdel(pcols,pver):pressure thickness of levels (Pa)
  // zm(pcols,pver)  :altitude (above ground) at level centers (m)
  // pblh(pcols)     :planetary boundary layer depth (m)
  // qv(pcols,pver)  :specific humidity (kg/kg)
  // cld(ncol,pver)  :cloud fraction (-) *** NOTE ncol dimension
  // dgncur_a(pcols,pver,ntot_amode)
  // dgncur_awet(pcols,pver,ntot_amode)
  //                                :dry & wet geo. mean dia. (m) of
  // number distrib. wetdens_host(pcols,pver,ntot_amode): interstitial
  // aerosol wet density (kg/m3)

  // DESCRIPTION:
  // calculates changes to gas and aerosol TMRs (tracer mixing ratios) from
  //    gas-aerosol exchange (condensation/evaporation)
  //    growth from smaller to larger modes (renaming) due to both
  //       condensation and cloud chemistry
  //    new particle nucleation
  //    coagulation
  //    transfer of particles from hydrophobic modes to hydrophilic modes
  //    (aging)
  //       due to condensation and coagulation
  //
  // the incoming mixing ratios (qq and qqcw) are updated before output
  //
  // REVISION HISTORY:
  //   RCE 07.04.13:  Adapted from earlier version of CAM5 modal aerosol
  //   routines
  //                  for these processes
  //

  // qgcmN and qqcwgcmN (N=1:4) are grid-cell mean tracer mixing ratios
  // (TMRs, mol/mol or #/kmol)
  //    N=1 - before gas-phase chemistry
  //    N=2 - before cloud chemistry
  //    N=3 - incoming values (before gas-aerosol exchange, newnuc, coag)
  //    N=4 - outgoing values (after  gas-aerosol exchange, newnuc, coag)

  // qsubN and qqcwsubN (N=1:4) are TMRs in sub-areas
  //    currently there are just clear and cloudy sub-areas
  //    the N=1:4 have same meanings as for qgcmN

  // Compute saturation vapor pressure
  const Real epsqs = haero::Constants::weight_ratio_h2o_air;

  // Saturation vapor pressure
  const Real ev_sat = conversions::vapor_saturation_pressure_magnus(temp, pmid);

  // Saturation specific humidity
  const Real qv_sat = epsqs * ev_sat / (pmid - (1 - epsqs) * ev_sat);

  // total # of subareas to do calculations for
  int nsubarea;

  // total # of cloudy subareas
  int ncldy_subarea;

  // indices of the clear and cloudy subareas
  int jclea, jcldy;

  // whether a subarea is cloudy
  bool iscldy_subarea[maxsubarea()];

  // area fraction of each active subarea[unitless]
  Real afracsub[maxsubarea()];

  // cloudy and clear fractions of the grid cell
  Real fcldy, fclea;

  setup_subareas(cld,                                      // in
                 nsubarea, ncldy_subarea, jclea, jcldy,    // out
                 iscldy_subarea, afracsub, fclea, fcldy);  // out
  EKAT_KERNEL_ASSERT_MSG(nsubarea < maxsubarea(),
                         "Error! modal_aero_amicphys_intr: "
                         "nsubarea should be < maxsubarea() \n");

  const Real relhumgcm = haero::max(0.0, haero::min(1.0, qv / qv_sat));

  Real relhumsub[maxsubarea()];
  set_subarea_rh(ncldy_subarea, jclea, jcldy, afracsub, relhumgcm,  // in
                 relhumsub);                                        // out

  //-------------------------------
  // Set aerosol water in subareas
  //-------------------------------
  // Notes from Dick Easter/Steve Ghan: how to treat aerosol water in
  // subareas needs more work/thinking Currently modal_aero_water_uptake
  // calculates qaerwat using the grid-cell mean interstital-aerosol
  // mix-rates and the clear-area RH. aerosol water mixing ratios (mol/mol)
  Real qaerwatsub3[num_modes][maxsubarea()];
  assign_2d_array(num_modes, maxsubarea(), 0.0,  // in
                  qaerwatsub3);                  // out

  //-------------------------------------------------------------------------
  // Set gases, interstitial aerosols, and cloud-borne aerosols in subareas
  //-------------------------------------------------------------------------
  // Copy grid cell mean mixing ratios; clip negative values if any.
  Real qgcm1[gas_pcnst], qgcm2[gas_pcnst], qgcm3[gas_pcnst];
  Real qqcwgcm2[gas_pcnst], qqcwgcm3[gas_pcnst];  // cld borne aerosols
  for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
    // Gases and interstitial aerosols
    qgcm1[icnst] = haero::max(0, q_pregaschem[icnst]);
    qgcm2[icnst] = haero::max(0, q_precldchem[icnst]);
    qgcm3[icnst] = haero::max(0, qq[icnst]);

    // Cloud-borne aerosols
    qqcwgcm2[icnst] = haero::max(0, qqcw_precldchem[icnst]);
    qqcwgcm3[icnst] = haero::max(0, qqcw[icnst]);
  }
  // Partition grid cell mean to subareas
  Real qsub1[gas_pcnst][maxsubarea()];
  Real qsub2[gas_pcnst][maxsubarea()];
  Real qsub3[gas_pcnst][maxsubarea()];
  Real qqcwsub1[gas_pcnst][maxsubarea()];
  Real qqcwsub2[gas_pcnst][maxsubarea()];
  Real qqcwsub3[gas_pcnst][maxsubarea()];

  set_subarea_gases_and_aerosols(nsubarea, jclea, jcldy, fclea,  // in
                                 fcldy,                          // in
                                 qgcm1, qgcm2, qqcwgcm2, qgcm3,  // in
                                 qqcwgcm3,                       // in
                                 qsub1, qsub2, qqcwsub2, qsub3,  // out
                                 qqcwsub3);                      // out

  //  Initialize the "after-amicphys" values
  Real qsub4[gas_pcnst][maxsubarea()]       = {};
  Real qqcwsub4[gas_pcnst][maxsubarea()]    = {};
  Real qaerwatsub4[num_modes][maxsubarea()] = {};

  //
  // start integration
  //
  Real dgn_a[num_modes]    = {};
  Real dgn_awet[num_modes] = {};
  Real wetdens[num_modes]  = {};
  assign_1d_array(num_modes, 0.0,     // in
                  dgn_a);             // out
  assign_1d_array(num_modes, 0.0,     // in
                  dgn_awet);          // out
  assign_1d_array(num_modes, 1000.0,  // in
                  wetdens);           // out

  for(int imode = 0; imode < num_modes; ++imode) {
    dgn_a[imode]    = dgncur_a[imode];
    dgn_awet[imode] = dgncur_awet[imode];
    wetdens[imode]  = haero::max(1000.0, wetdens_host[imode]);
  }

  Real qsub_tendaa[gas_pcnst][nqtendaa()][maxsubarea()]       = {};
  Real qqcwsub_tendaa[gas_pcnst][nqqcwtendaa()][maxsubarea()] = {};
  mam_amicphys_1gridcell(
      // in
      config, deltat, nsubarea, ncldy_subarea, iscldy_subarea, afracsub, temp,
      pmid, pdel, zm, pblh, relhumsub, dgn_a, dgn_awet, wetdens, qsub1, qsub2,
      qqcwsub2, qsub3, qqcwsub3,
      // inout
      qaerwatsub3, qsub4, qqcwsub4, qaerwatsub4, qsub_tendaa, qqcwsub_tendaa);

  //=================================================================================================
  // Aerosol microphysics calculations done for all subareas. Form new grid cell
  // mean mixing ratios.
  //=================================================================================================
  // Gases and aerosols
  //----------------------
  // Calculate new grid cell mean values
  Real qgcm4[gas_pcnst];
  Real qqcwgcm4[gas_pcnst];

  form_gcm_of_gases_and_aerosols_from_subareas(
      // in
      nsubarea, ncldy_subarea, afracsub, qsub4, qqcwsub4, qqcwgcm3,
      // out
      qgcm4, qqcwgcm4);

  // Copy grid cell mean values to output arrays
  for(int icnst = 0; icnst < gas_pcnst; ++icnst) {
    if(lmapcc_all(icnst) > 0) {
      qq[icnst] = qgcm4[icnst];
    }
    if(lmapcc_all(icnst) >= lmapcc_val_aer()) {
      qqcw[icnst] = qqcwgcm4[icnst];
    }
  }

  //================================================================
  // Process diagnostics of the current grid cell
  //================================================================
  Real qgcm_tendaa[gas_pcnst][nqtendaa()];
  Real qqcwgcm_tendaa[gas_pcnst][nqqcwtendaa()];
  get_gcm_tend_diags_from_subareas(
      // in
      nsubarea, ncldy_subarea, afracsub, qsub_tendaa, qqcwsub_tendaa,
      // out
      qgcm_tendaa, qqcwgcm_tendaa);

#if 0
  //This code is for diagnostics only
  // Get gravity
  using C                      = physics::Constants<Real>;
  static constexpr auto gravit = C::gravit;  // Gravity [m/s2]

accumulate_column_tend_integrals( pdel, gravit,                         // in
                                        qgcm_tendaa,         qqcwgcm_tendaa,         // in
                                        q_coltendaa(ii,:,:), qqcw_coltendaa(ii,:,:)  )// inout

ncluster_3dtend_nnuc(ii,kk) = misc_vars_aa%ncluster_tend_nnuc_1grid

#endif
}  // modal_aero_amicphys_intr

}  // namespace scream::impl
