#ifndef SCREAM_UTILS_HPP
#define SCREAM_UTILS_HPP

#include "share/scream_types.hpp"

#include <ekat/ekat_assert.hpp>
#include <ekat/kokkos/ekat_kokkos_types.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <iterator>
#include <list>
#include <algorithm>
#include <map>
#include <iostream>

namespace scream {

enum MemoryUnits {
  B = 1,
  KB,
  MB,
  GB,
  KiB,
  MiB,
  GiB
};

template<typename VT>
typename VT::HostMirror
cmvdc (const VT& v)
{
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

// Gets current memory (RAM) usage by current process.
long long get_mem_usage (const MemoryUnits u);

// Micro-utility, that given an enum returns the underlying int.
// The only use of this is if you need to sort scoped enums.
template<typename EnumT>
KOKKOS_FUNCTION
constexpr typename
std::enable_if<std::is_enum<EnumT>::value,
               typename std::underlying_type<EnumT>::type
              >::type
etoi (const EnumT e) {
  return static_cast<typename std::underlying_type<EnumT>::type>(e);
}

inline void broadcast_string (std::string& s, const ekat::Comm& comm, const int root) {
  int size = s.size();
  comm.broadcast(&size,1,root);
  s.resize(size);
  comm.broadcast(&s.front(),size,root);
}

// Utility function, to work around a funcky gcc8+cuda10 issue,
// where calling sort() on a length-2 list returns a length-1 list.
template<typename T>
void sort (std::list<T>& l) {
  if (l.size()==2) {
    if (l.back()<l.front()) {
      std::swap(l.front(),l.back());
    }
    return;
  }
  l.sort();
}

// This routine tries to find an arrangment of elements that allow each
// of the input groups to be a contiguous subarray of the global arrangement.
// E.g., given the groups of elments
//    G1=(A,B,C), G2=(A,B,C,D,E), G3=(C,D), G4=(C,D,E,F), G5=((D,E,F,G).
// the ordering of elements (A,B,C,D,E,F,G) clearly has all G's as contiguous
// subset, though it might not be obvious if one scrambled the entris within
// each input group.
// If such an ordering cannot be found, we return an empty list.

// The algorithm is the following.
//  - We proceed iteratively, trying to add one group at a time to the result.
//  - Internally, we store the result as a list of lists (LoL), like
//        [ [A,B], [C], [D,E,F], [G,H], [I], ... ]
//    If the algorithm succeeds, the output of the function is a single list
//    obtainied by splicing together all the inner lists.
//  - We are only allowed to modify the outer list in the following ways:
//      a) rearrange the elements in an inner list: [D,E,F]->[F,D,E]
//      b) split an inner list into multiple lists: [[A,B,C],[D,E]]->[[A,B],[C],[D,E}]
//      c) append a new list (that does not overlap with existing sublists):
//        [[A,B],[C,D]] -> [[A,B],[C,D],[E,F,G]]
//  - The 1st group is added as a single sublist, that is, given G1,..,G5 above,
//    after the first iteration, the list of lists would be [[A,B,C]].
//  - At the generic k-th step of the algorithm, we try to add group G to the LoL,
//    doing the following:
//     - compute intersection of G with each of the sublists of LoL
//     - check that all non-empty intersections are 'contiguous'. If not, the
//       overall algorithm failed, and we can return.
//     - if there's a part of G that does not intersect with any sublist, then
//       there are two scenarios:
//        i) all intersections are at the beginning or end of LoL: good.
//        ii) intersections are neither at the beginning nor at the end of LoL.
//            In this case, the overall algorithm failed, and we can return.
//     - if G has a non-empty intersection with sublist S, there are two scenarios:
//        i) the intersection is the whole S: nothing to do.
//        ii) the intersection is smaller than S: two sub-scenarios:
//          *) S is the first or last sublist with non-empty intersection:
//             split S into [intersection][rest] or [rest][intersection]
//             (the former if it's the first intersection, the former if it's the last).
//             This would correspond to operation (a+b) on the LoL.
//          *) S is neither the first nor the last. In this case the overall
//             algorithm has failed, and we can return.
//     - if G has a portion that does not intersect with any sublist S, add it
//       as a new sublist at the end. This corresponds to operation (c) from above.

// Regardless of the order in which we process the individual lists, if there is
// an arrangement A of T's that allows to have G1,...,Gn as a contiguous sublist
// of A, the above algorithm is guaranteed to find it.
template<typename T>
std::list<T> contiguous_superset (const std::list<std::list<T>>& groups)
{
  // Avoid accessing empty lists.
  if (groups.size()==0) {
    return std::list<T>();
  }

  // Intersection and difference of two lists. Simply wrap std::set_intersection
  // and std::set_difference resepectively, allowing a lighter syntax.
  auto intersect = [] (const std::list<T>& lhs,
                       const std::list<T>& rhs) -> std::list<T> {
    std::list<T> out;
    std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                          std::back_inserter(out));
    return out;
  };
  auto difference = [] (const std::list<T>& lhs,
                       const std::list<T>& rhs) -> std::list<T> {
    std::list<T> out;
    std::set_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                        std::back_inserter(out));
    return out;
  };

  // Would be nice if list::sort and list::unique returned the list,
  // so one could chain them into operations. Since they don't, use lambdas.
  // Also, we don't want to modify input list just to check uniqueness,
  // so create copies.
  auto sort = [] (const std::list<T>& l) -> std::list<T> {
    auto copy = l;
    ::scream::sort(copy);
    return copy;
  };
  auto unique = [] (const std::list<T>& l) -> std::list<T> {
    auto copy = l;
    copy.unique();
    return copy;
  };

  // Check inputs: we require unique lists, already sorted.
  for (const auto& g : groups) {
    EKAT_REQUIRE_MSG(sort(g)==g,
        "Error! Individual input lists must already be sorted.\n");
    EKAT_REQUIRE_MSG(unique(sort(g)).size()==g.size(),
        "Error! Individual input lists must not contain repeated elements.\n");
  }

  // A list-of-lists (lol);
  std::list<std::list<T>> lol;

  decltype(lol.begin()) it_lol;
  for (const auto& g : groups) {
    // Intersect this group with each inner list inside lol, and keep track of
    // where non-empty intersections happen. These intersections must happen
    // with a contiguous set of inner lists (othewise the algo fails, see
    // description above). Store interval where they happen as [first_pos,last_pos) range.
    // Note: \cap is the tex command for the set intersection symbol.
    std::vector<std::list<T>> caps;
    std::vector<size_t> caps_pos;
    auto remainder = g;
    for (auto it=lol.begin(); it!=lol.end(); ++it) {
      auto tmp = intersect(*it,g);
      if (tmp.size()>0) {
        // Found non-empty intersection. Track it.
        caps.emplace_back(std::move(tmp));
        remainder = difference(remainder,caps.back());
        caps_pos.push_back(std::distance(lol.begin(),it));
      }
    }

    // If we don't have any intersection, append this group to lol, and continue
    // Note: this for sure happens when processing the first group
    if (caps_pos.size()==0) {
      lol.push_back(g);
      continue;
    }

    // If we have caps, then they must be with a contiguous set of inner list
    // (otherwise this group would be fragmented).
    // This can be checked easily by inspecting 1st and last entry in caps_pos
    if ( (caps_pos.back()-caps_pos.front()+1)!=caps_pos.size()) {
      return std::list<T>();
    }

    // If we have a non-empty reminder, then all the caps must be
    // either at the front or the back of the lol, that is, either the first intersection
    // is at pos=0, or the last at pos=lol.size()-1, or both (otherwise this group would be fragmented).
    if (remainder.size()>0 && !(caps_pos.front()==0 || caps_pos.back()==(lol.size()-1))) {
      return std::list<T>();
    }

    // If there is a remainder, then either the first or the last intersection must
    // be "complete", meaning that the intersection is the whole inner list.
    // If not, there would be fragmentation. E.g., if G=[B,D,E], and
    // lol=[[A,B],[C,D]], there's no rearrangement that works.
    if (remainder.size()>0 && caps_pos.size()>1) {
      auto it_lol_first = std::next(lol.begin(),caps_pos.front());
      auto it_lol_last  = std::next(lol.begin(),caps_pos.back());
      if (it_lol_first->size()>caps.front().size() &&
          it_lol_last->size()>caps.back().size()) {
        return std::list<T>();
      }
    }

    // If we have 3+ caps, only the first and last caps can be
    // less than the corresponding inner list in lol.
    for (size_t i=1; i<(caps_pos.size()-1); ++i) {
      it_lol = std::next(lol.begin(),caps_pos[i]);
      const auto& intersection = caps[i];
      if (it_lol->size()!=intersection.size()) {
        return std::list<T>();
      }
    }

    // Ok, we successfully passed all checks. Now we have to do two things:
    //  - insert the remainder (if any) at the front or back
    //  - if first/last caps are smaller than corresponding inner list,
    //    split the inner list in two
    it_lol = std::next(lol.begin(),caps_pos.front());
    if (caps.front().size()<it_lol->size()) {
      // The first intersection is smaller than the inner list.
      // Split the inner list in intersection and diff=list-intersection.
      // The order in which the two sublists need to be added
      // depends on whether there is any other intersection after this:
      //  - more caps: [diff],[intersection]
      //  - no more caps: [intersection],[diff]
      // These choices allow to keep the current group contiguous.
      // Note: the latter is really needed only if remainder is non-empty.

      auto diff = difference(*it_lol,caps.front());
      if (caps_pos.size()>1) {
        *it_lol = std::move(caps.front());
        lol.emplace(it_lol,std::move(diff));
      } else {
        *it_lol = std::move(diff);
        lol.emplace(it_lol,std::move(caps.front()));
      }

      // Note: since we added something right before where the 1st
      //       intersection happened, all caps_pos indices need to
      //       be updated (adding 1).
      for (size_t i=1; i<caps_pos.size(); ++i) {
        ++caps_pos[i];
      }
    }

    // Note: process caps_pos.back() only if there's 2+ caps
    //       (otherwise we process the same intersection twice)
    it_lol = std::next(lol.begin(),caps_pos.back());
    if (caps_pos.size()>1 && caps.back().size()<it_lol->size()) {
      // The last intersection is smaller than the inner list.
      // Split the inner list in intersection and diff=list-intersection.
      // Unlike the previous case, here we know that there were 2+ caps,
      // So we can add the two lists as [intersection][diff].
      // This will allow the current group to remain contiguous.
      auto diff = difference(*it_lol,caps.back());
      *it_lol = std::move(diff);
      lol.emplace(it_lol,std::move(caps.back()));
    }


    if (remainder.size()>0) {
      if (caps_pos.front()==0) {
        lol.push_front(remainder);
      } else {
        lol.push_back(remainder);
      }
    }
  }

  // We made it! Now splice all sublists together.
  std::list<T> out;
  for (auto& l : lol) {
    out.splice(out.end(),l);
  }
  return out;
}

/* Given a column of data for variable "label" from the reference run
 * (probably master) and from your new exploratory run, loop over all
 * heights and confirm whether or not the relative difference between
 * runs is within tolerance "tol". If not, print debug info. Here, "a"
 * is the value from the reference run and "b" is from the new run.
 * This is used by the run_and_cmp tests.
 */
template <typename Scalar, typename Toltype>
Int compare (const std::string& label, const Scalar* a,
             const Scalar* b, const Int& n, const Toltype& tol) {

  Int nerr1 = 0;
  Int nerr2 = 0;
  Scalar den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Scalar worst = 0;
  for (Int i = 0; i < n; ++i) {
    if (std::isnan(a[i]) || std::isinf(a[i]) ||
        std::isnan(b[i]) || std::isinf(b[i])) {
      ++nerr1;
      continue;
    }

    const auto num = std::abs(a[i] - b[i]);
    if (num > tol*den) {
      ++nerr2;
      worst = std::max(worst, num);
    }
  }

  if (nerr1) {
    std::cout << label << " has " << nerr1 << " infs + nans.\n";

  }

  if (nerr2) {
    std::cout << label << " > tol " << nerr2 << " times. Max rel diff= " << (worst/den)
             << " normalized by ref impl val=" << den << ".\n";

  }

  return nerr1 + nerr2;
}

inline void
check_mpi_call (int err, const std::string& context) {
  EKAT_REQUIRE_MSG (err==MPI_SUCCESS,
      "Error! MPI operation encountered an error.\n"
      "  - err code: " + std::to_string(err) + "\n"
      "  - context: " + context + "\n");
}

// Find the full filename list from patterns
std::vector<std::string> filename_glob(const std::vector<std::string>& patterns);

// Use globloc for each filename pattern
std::vector<std::string> globloc(const std::string& pattern);

constexpr int eamxx_swbands() {
  // This function returns the total number of SW bands in RRTMGP,
  return 14;
}

constexpr int eamxx_vis_swband_idx() {
  // This function returns the index of the visible SW band in RRTMGP,
  // which currently (as of 2024-04-23) is supposed to be 10.
  // This index is used in the AODVis diagnostic, and should ideally
  // be shared across interested processes for further diagnostics.
  // This index (10) corresponds to the band that has wavelength 550 nm.
  return 10;
}

struct DefaultMetadata {

  std::string get_longname (const std::string& name) {
    if (name_2_longname.count(name)>0) {
      return name_2_longname.at(name);
    } else {
      // TODO: Do we want to print a Warning message?  I'm not sure if its needed.
      return name;
    }
  }

  std::string get_standardname (const std::string& name) {
    if (name_2_standardname.count(name)>0) {
      return name_2_standardname.at(name);
    } else {
      // TODO: Do we want to print a Warning message?  I'm not sure if its needed.
      return name;
    }
  }

  // Create map of longnames, can be added to as developers see fit.
  std::map<std::string,std::string> name_2_longname = {
    {"lev","hybrid level at midpoints (1000*(A+B))"},
    {"ilev","hybrid level at interfaces (1000*(A+B))"},
    {"hyai","hybrid A coefficient at layer interfaces"},
    {"hybi","hybrid B coefficient at layer interfaces"},
    {"hyam","hybrid A coefficient at layer midpoints"},
    {"hybm","hybrid B coefficient at layer midpoints"}
  };

  // Create map of longnames, can be added to as developers see fit.
  std::map<std::string,std::string> name_2_standardname = {
    {"p_mid"                                                       , "air_pressure"},
    {"p_mid_at_cldtop"                                             , "air_pressure_at_cloud_top"},
    {"T_2m"                                                        , "air_temperature"},
    {"T_mid"                                                       , "air_temperature"},
    {"T_mid_at_cldtop"                                             , "air_temperature_at_cloud_top"},
    {"aero_g_sw"                                                   , "asymmetry_factor_of_ambient_aerosol_particles"},
    {"pbl_height"                                                  , "atmosphere_boundary_layer_thickness"},
    {"precip_liq_surf_mass"                                        , "atmosphere_mass_content_of_liquid_precipitation"},
    {"cldlow"                                                      , "low_type_cloud_area_fraction"},
    {"cldmed"                                                      , "medium_type_cloud_area_fraction"},
    {"cldhgh"                                                      , "high_type_cloud_area_fraction"},
    {"cldtot"                                                      , "cloud_area_fraction"},
    {"cldfrac_tot_at_cldtop"                                       , "cloud_area_fraction"},
    {"cldfrac_tot"                                                 , "cloud_area_fraction_in_atmosphere_layer"},
    {"cldfrac_tot_for_analysis"                                    , "cloud_area_fraction_in_atmosphere_layer"},
    {"cldfrac_rad"                                                 , "cloud_area_fraction_in_atmosphere_layer"},
    {"qi"                                                          , "cloud_ice_mixing_ratio"},
    {"qc"                                                          , "cloud_liquid_water_mixing_ratio"},
    {"U"                                                           , "eastward_wind"},
    {"eff_radius_qi"                                               , "effective_radius_of_cloud_ice_particles"},
    {"eff_radius_qc"                                               , "effective_radius_of_cloud_liquid_water_particles"},
    {"eff_radius_qc_at_cldtop"                                     , "effective_radius_of_cloud_liquid_water_particles_at_liquid_water_cloud_top"},
    {"eff_radius_qr"                                               , "effective_radius_of_cloud_rain_particles"},
    {"qv"                                                          , "humidity_mixing_ratio"},
    {"cldfrac_ice_at_cldtop"                                       , "ice_cloud_area_fraction"},
    {"cldfrac_ice"                                                 , "ice_cloud_area_fraction_in_atmosphere_layer"},
    {"omega"                                                       , "lagrangian_tendency_of_air_pressure"},
    {"landfrac"                                                    , "land_area_fraction"},
    {"latitude"                                                    , "latitude"},
    {"cldfrac_liq_at_cldtop"                                       , "liquid_water_cloud_area_fraction"},
    {"cldfrac_liq"                                                 , "liquid_water_cloud_area_fraction_in_atmosphere_layer"},
    {"longitude"                                                   , "longitude"},
    {"rainfrac"                                                    , "mass_fraction_of_liquid_precipitation_in_air"},
    {"V"                                                           , "northward_wind"},
    {"nc"                                                          , "number_concentration_of_cloud_liquid_water_particles_in_air"},
    {"cdnc_at_cldtop"                                              , "number_concentration_of_cloud_liquid_water_particles_in_air_at_liquid_water_cloud_top"},
    {"ni"                                                          , "number_concentration_of_ice_crystals_in_air"},
    {"aero_tau_sw"                                                 , "optical_thickness_of_atmosphere_layer_due_to_ambient_aerosol_particles"},
    {"aero_tau_lw"                                                 , "optical_thickness_of_atmosphere_layer_due_to_ambient_aerosol_particles"},
    {"aero_ssa_sw"                                                 , "single_scattering_albedo_in_air_due_to_ambient_aerosol_particles"},
    {"sunlit"                                                      , "sunlit_binary_mask"},
    {"ps"                                                          , "surface_air_pressure"},
    {"LW_flux_dn_at_model_bot"                                     , "surface_downwelling_longwave_flux_in_air"},
    {"SW_flux_dn_at_model_bot"                                     , "surface_downwelling_shortwave_flux_in_air"},
    {"SW_clrsky_flux_dn_at_model_bot"                              , "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky"},
    {"phis"                                                        , "surface_geopotential"},
    {"surf_radiative_T"                                            , "surface_temperature"},
    {"surf_sens_flux"                                              , "surface_upward_sensible_heat_flux"},
    {"SW_flux_dn_at_model_top"                                     , "toa_incoming_shortwave_flux"},
    {"LW_flux_up_at_model_top"                                     , "toa_outgoing_longwave_flux"},
    {"LW_clrsky_flux_up_at_model_top"                              , "toa_outgoing_longwave_flux_assuming_clear_sky"},
    {"surf_evap"                                                   , "water_evapotranspiration_flux"},
    {"AtmosphereDensity"                                           , "air_density"},
    {"PotentialTemperature"                                        , "air_potential_temperature"},
    {"SeaLevelPressure"                                            , "air_pressure_at_mean_sea_level"},
    {"IceWaterPath"                                                , "atmosphere_mass_content_of_cloud_ice"},
    {"LiqWaterPath"                                                , "atmosphere_mass_content_of_cloud_liquid_water"},
    {"VapWaterPath"                                                , "atmosphere_mass_content_of_water_vapor"},
    {"AerosolOpticalDepth550nm"                                    , "atmosphere_optical_thickness_due_to_ambient_aerosol_particles"},
    {"Exner"                                                       , "dimensionless_exner_function"},
    {"z_mid"                                                       , "geopotential_height"},
    {"geopotential_mid"                                            , "geopotential_height"},
    {"RelativeHumidity"                                            , "relative_humidity"},
    {"surface_upward_latent_heat_flux"                             , "surface_upward_latent_heat_flux"},
    {"LongwaveCloudForcing"                                        , "toa_longwave_cloud_radiative_effect"},
    {"ShortwaveCloudForcing"                                       , "toa_shortwave_cloud_radiative_effect"},
    {"VirtualTemperature"                                          , "virtual_temperature"},
    {"VaporFlux"                                                   , "water_evapotranspiration_flux"},
    {"wind_speed"                                                  , "wind_speed"}
  };
  
};


} // namespace scream

#endif // SCREAM_UTILS_HPP
