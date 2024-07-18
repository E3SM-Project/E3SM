#ifndef MAM_EMISSIONS_READ_TABLES_HPP
#define MAM_EMISSIONS_READ_TABLES_HPP

#include "ekat/ekat_parameter_list.hpp"
#include "mam_coupling.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

// later to mam_coupling.hpp
namespace scream::mam_coupling {

using view_1d_host = typename KT::view_1d<Real>::HostMirror;
using view_1d_int_host = typename KT::view_1d<int>::HostMirror;
using view_2d_host = typename KT::view_2d<Real>::HostMirror;
// using view_5d_host    = typename KT::view_ND<Real,5>::HostMirror;
// using complex_view_1d = typename KT::view_1d<Kokkos::complex<Real>>;

constexpr int n_srf_emiss = mam4::mo_srf_emissions::n_srf_emiss;
constexpr int n_online_emiss = mam4::aero_model_emissions::n_online_emiss;

// FIXME: this will need to change when we remap to flattened column idx
constexpr int nlat_srf = 96;
constexpr int	nlon_srf = 144;

constexpr int nalti_online = 13;
constexpr int nlat_online = 96;
constexpr int nlon_online = 144;

using namespace ShortFieldTagsNames;

// std::map<std::string, std::vector<std::string>> map_srf_emiss_file_vars;

inline void set_file_var_names(std::map<std::string, std::vector<std::string>> &var_map,
                               std::map<std::string, int> &spec_map) {
  // for (const auto &spec : spec_map) {
  //   std::string spec_name = spec.first;
  //   std::cout << "var_map[spec_name] = " << var_map[spec_name] << "\n";
  // }
}

// struct AerosolSurfaceEmissionsHostData {
//   // these have dim = n_species
//   view_1d_host emis_species_index;
//   view_1d_host emis_species_units;
//   view_1d_host emis_species_name;
//   // molecular weight
//   view_1d_host emis_species_mw;
//   // number of sectors in each field
//   view_1d_int_host emis_species_nsectors;
//   // FIXME: not quite sure what this does--maybe just a placeholder for
//   // fields(:, i_sector)?
//   view_1d_host emis_species_sector;
//   // note fields have dim = n_species x nsectors
//   // TODO: fields have units??? maybe the same as the upper spec units
//   view_2d_host emis_species_fields;
// };

// using AerosolSurfaceEmissionsDeviceData =
    // mam4::mo_srf_emissions::AerosolSurfaceEmissionsDeviceData;

// inline void set_emissions_params(
//     AerosolSurfaceEmissionsHostData& aerosol_emissions_host_data,
//     ekat::ParameterList& params_emissions,
//     std::map<std::string, FieldLayout>& layouts,
//     std::map<std::string, view_1d_host>& host_views) {
//   // Set up input structure to read data from file.
//   using strvec_t = std::vector<std::string>;
//   // using namespace ShortFieldTagsNames;

//   // using SrfEmisDims =
//   // mam4::mo_srf_emissions::AerosolSurfaceEmissionsDimensions; SrfEmisDims
//   // srf_emimssions_dims;
// }

// inline void set_emissions_names(const std::map<std::string, int> map_spec_id,
//                     const std::string emis_type,
//                     const ekat::ParameterList& m_params,
//                     std::map<std::string, view_1d_host>& host_views) {

//   using view_1d_host = typename KT::view_1d<Real>::HostMirror;

//   std::string

// } // end set_emissions_names

} // namespace scream::mam_coupling

#endif
