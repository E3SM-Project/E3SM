#include <physics/mam/eamxx_mam_generic_process_interface.hpp>

namespace scream {
// ================================================================
//  Constructor
// ================================================================

MAMGenericInterface::MAMGenericInterface(const ekat::Comm &comm,
                                           const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}
// ================================================================
void MAMGenericInterface::add_invariant_check_for_aerosol(
  //std::shared_ptr<const AbstractGrid> grid_
  )
{
// NOTE: Using only one range for all num variables.
  const std::string nmr_label="nmr";
  const auto min_nmr = mam_coupling::physical_min(nmr_label);
  const auto max_nmr = mam_coupling::physical_max(nmr_label);
  const std::string mmr_label="mmr";
  const auto min_mmr = mam_coupling::physical_min(mmr_label);
  const auto max_mmr = mam_coupling::physical_max(mmr_label);

  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    const std::string int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(mode);

    // NOTE: invariant_check checks both add_precondition_check and add_postcondition_check
    // Because these variables are updated, we will check them before and after computations.
    add_invariant_check<FieldWithinIntervalCheck>(
    get_field_out(int_nmr_field_name),grid_,min_nmr,max_nmr,false);

   const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    add_invariant_check<FieldWithinIntervalCheck>(
     get_field_out(cld_nmr_field_name),grid_,min_nmr,max_nmr,false);

   //
   for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios

      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(not int_mmr_field_name.empty()) {
         add_invariant_check<FieldWithinIntervalCheck>(
         get_field_out(int_mmr_field_name),grid_,min_mmr,max_mmr,false);
      }
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
      // NOT advected
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(mode, a);
      if(not cld_mmr_field_name.empty()) {
        add_invariant_check<FieldWithinIntervalCheck>(
        get_field_out(cld_mmr_field_name),grid_,min_mmr,max_mmr,false);
      }
    }  // end for loop num species
  }

  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_invariant_check<FieldWithinIntervalCheck>(
    get_field_out(gas_mmr_field_name),grid_,min_mmr,max_mmr,false);
   }  // end for loop num gases
}
// ================================================================
}  // namespace scream
