#ifndef MARINE_ORGANICS_IMPL_HPP
#define MARINE_ORGANICS_IMPL_HPP

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_timing.hpp"

namespace scream {
namespace marine_organics {

template <typename S, typename D>
std::shared_ptr<AbstractRemapper>
marineOrganicsFunctions<S, D>::create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &data_file, const std::string &map_file,
    const std::vector<std::string> &field_name, const std::string &dim_name1) {
  using namespace ShortFieldTagsNames;

  scorpio::register_file(data_file, scorpio::Read);
  const int ncols_data = scorpio::get_dimlen(data_file, dim_name1);

  scorpio::release_file(data_file);

  // Since shallow clones are cheap, we may as well do it (less lines of
  //  code)
  auto horiz_interp_tgt_grid =
      model_grid->clone("marine_organics_horiz_interp_tgt_grid", true);

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if(ncols_data == ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(ncols_data <= ncols_model,
                     "Error! We do not allow to coarsen marine organics "
                     "data to fit the model. We only allow\n"
                     "       marine organics data to be at the same or "
                     "coarser resolution as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(map_file != "",
                     "ERROR: marine organics data is on a different grid "
                     "than the model one,\n"
                     "       but remap file is missing from marine organics "
                     "parameter list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();

  const auto tgt_grid = remapper->get_tgt_grid();

  const auto layout_2d = tgt_grid->get_2d_scalar_layout();
  using namespace ekat::units;
  using namespace ekat::prefixes;
  Units umolC(micro * mol, "umol C");

  std::vector<Field> fields_vector;

  const int field_size = field_name.size();
  for(int icomp = 0; icomp < field_size; ++icomp) {
    auto comp_name = field_name[icomp];
    // set and allocate fields
    Field f(FieldIdentifier(comp_name, layout_2d, umolC, tgt_grid->name()));
    f.allocate_view();
    fields_vector.push_back(f);
    remapper->register_field_from_tgt(f);
  }

  remapper->registration_ends();

  return remapper;

}  // create_horiz_remapper

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
std::shared_ptr<AtmosphereInput>
marineOrganicsFunctions<S, D>::create_data_reader(
    const std::shared_ptr<AbstractRemapper> &horiz_remapper,
    const std::string &data_file) {
  std::vector<Field> io_fields;
  for(int ifld = 0; ifld < horiz_remapper->get_num_fields(); ++ifld) {
    io_fields.push_back(horiz_remapper->get_src_field(ifld));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(data_file, io_grid, io_fields, true);
}  // create_data_reader

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
void marineOrganicsFunctions<S, D>::update_marine_organics_data_from_file(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    const int &time_index,  // zero-based
    AbstractRemapper &horiz_interp, marineOrganicsInput &marineOrganics_input) {
  start_timer("EAMxx::marineOrganics::update_marine_organics_data_from_file");

  // 1. Read from file
  start_timer(
      "EAMxx::marineOrganics::update_marine_organics_data_from_file::read_"
      "data");
  scorpio_reader->read_variables();
  stop_timer(
      "EAMxx::marineOrganics::update_marine_organics_data_from_file::read_"
      "data");

  // 2. Run the horiz remapper (it is a do-nothing op if marineOrganics data is
  // on same grid as model)
  start_timer(
      "EAMxx::marineOrganics::update_marine_organics_data_from_file::horiz_"
      "remap");
  horiz_interp.remap(/*forward = */ true);
  stop_timer(
      "EAMxx::marineOrganics::update_marine_organics_data_from_file::horiz_"
      "remap");

  // 3. Get the tgt field of the remapper
  start_timer(
      "EAMxx::marineOrganics::update_marine_organics_data_from_file::get_"
      "field");
  // Recall, the fields are registered in the order:
  // Read the field from the file

  for(int ifld = 0; ifld < horiz_interp.get_num_fields(); ++ifld) {
    auto sector = horiz_interp.get_tgt_field(ifld).get_view<const Real *>();
    const auto emiss = Kokkos::subview(marineOrganics_input.data.emiss_sectors,
                                       ifld, Kokkos::ALL());
    Kokkos::deep_copy(emiss, sector);
  }

  Kokkos::fence();

  stop_timer(
      "EAMxx::marineOrganics::update_marine_organics_data_from_file::get_"
      "field");

  stop_timer("EAMxx::marineOrganics::update_marine_organics_data_from_file");

}  // END update_marine_organics_data_from_file

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
void marineOrganicsFunctions<S, D>::update_marine_organics_timestate(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    AbstractRemapper &horiz_interp, marineOrganicsTimeState &time_state,
    marineOrganicsInput &beg, marineOrganicsInput &end) {
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that marineOrganics assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1;  // Make it 0-based
  if(month != time_state.current_month) {
    // Update the marineOrganics time state information
    time_state.current_month = month;
    time_state.t_beg_month =
        util::TimeStamp({ts.get_year(), month + 1, 1}, {0, 0, 0})
            .frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(), month + 1);

    // Copy end'data into beg'data, and read in the new
    // end
    std::swap(beg, end);

    // Update the marineOrganics forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong
    // values
    //       to be assigned.  A timestep greater than a month is very unlikely
    //       so we will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    update_marine_organics_data_from_file(scorpio_reader, ts, next_month,
                                          horiz_interp, end);
  }

}  // END updata_marine_organics_timestate

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
template <typename ScalarX, typename ScalarT>
KOKKOS_INLINE_FUNCTION ScalarX marineOrganicsFunctions<S, D>::linear_interp(
    const ScalarX &x0, const ScalarX &x1, const ScalarT &t) {
  return (1 - t) * x0 + t * x1;
}  // linear_interp

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
void marineOrganicsFunctions<S, D>::perform_time_interpolation(
    const marineOrganicsTimeState &time_state,
    const marineOrganicsInput &data_beg, const marineOrganicsInput &data_end,
    const marineOrganicsOutput &data_out) {
  using ExeSpace = typename KT::ExeSpace;
  using ESU      = ekat::ExeSpaceUtils<ExeSpace>;

  // Gather time stamp info
  auto &t_now   = time_state.t_now;
  auto &t_beg   = time_state.t_beg_month;
  auto &delta_t = time_state.days_this_month;

  // At this stage, begin/end must have the same dimensions
  EKAT_REQUIRE(data_end.data.ncols == data_beg.data.ncols);

  auto delta_t_fraction = (t_now - t_beg) / delta_t;

  EKAT_REQUIRE_MSG(delta_t_fraction >= 0 && delta_t_fraction <= 1,
                   "Error! Convex interpolation with coefficient out of "
                   "[0,1].\n  t_now  : " +
                       std::to_string(t_now) +
                       "\n"
                       "  t_beg  : " +
                       std::to_string(t_beg) +
                       "\n  delta_t: " + std::to_string(delta_t) + "\n");

  const int nsectors = data_beg.data.nsectors;
  const int ncols    = data_beg.data.ncols;
  using ExeSpace     = typename KT::ExeSpace;
  using ESU          = ekat::ExeSpaceUtils<ExeSpace>;
  const auto policy  = ESU::get_default_team_policy(ncols, nsectors);

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MemberType &team) {
        const int icol = team.league_rank();  // column index
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, 0u, nsectors), [&](int isec) {
              const auto beg = data_beg.data.emiss_sectors(isec, icol);
              const auto end = data_end.data.emiss_sectors(isec, icol);
              data_out.emiss_sectors(isec, icol) =
                  linear_interp(beg, end, delta_t_fraction);
            });
      });
  Kokkos::fence();

}  // perform_time_interpolation

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
void marineOrganicsFunctions<S, D>::marineOrganics_main(
    const marineOrganicsTimeState &time_state,
    const marineOrganicsInput &data_beg, const marineOrganicsInput &data_end,
    const marineOrganicsOutput &data_out) {
  // Beg/End/Tmp month must have all sizes matching

  EKAT_REQUIRE_MSG(
      data_end.data.ncols == data_beg.data.ncols,
      "Error! marineOrganicsInput data structs must have the same number of "
      "columns.\n");

  // Horiz interpolation can be expensive, and does not depend on the particular
  // time of the month, so it can be done ONCE per month, *outside*
  // marineOrganics_main (when updating the beg/end states, reading them from
  // file).
  EKAT_REQUIRE_MSG(data_end.data.ncols == data_out.ncols,
                   "Error! Horizontal interpolation is performed *before* "
                   "calling marineOrganics_main,\n"
                   "       marineOrganicsInput and marineOrganicsOutput data "
                   "structs must have the "
                   "same number columns "
                       << data_end.data.ncols << "  " << data_out.ncols
                       << ".\n");

  // Step 1. Perform time interpolation
  perform_time_interpolation(time_state, data_beg, data_end, data_out);
}  // marineOrganics_main

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
void marineOrganicsFunctions<S, D>::init_marine_organics_file_read(
    const int &ncol, const std::vector<std::string> &field_name,
    const std::string &dim_name1,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::string &data_file, const std::string &mapping_file,
    // output
    std::shared_ptr<AbstractRemapper> &marineOrganicsHorizInterp,
    marineOrganicsInput &data_start_, marineOrganicsInput &data_end_,
    marineOrganicsData &data_out_,
    std::shared_ptr<AtmosphereInput> &marineOrganicsDataReader) {
  // Init horizontal remap

  marineOrganicsHorizInterp = create_horiz_remapper(
      grid, data_file, mapping_file, field_name, dim_name1);

  // Initialize the size of start/end/out data structures
  data_start_ = marineOrganicsInput(ncol, field_name.size());
  data_end_   = marineOrganicsInput(ncol, field_name.size());
  data_out_.init(ncol, field_name.size(), true);

  // Create reader (an AtmosphereInput object)
  marineOrganicsDataReader =
      create_data_reader(marineOrganicsHorizInterp, data_file);

}  // init_marine_organics_file_read
}  // namespace marine_organics
}  // namespace scream

#endif  // MARINE_ORGANICS_IMPL_HPP