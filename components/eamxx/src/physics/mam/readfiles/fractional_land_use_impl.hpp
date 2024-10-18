#ifndef FRACTIONAL_LANDUSE_IMPL_HPP
#define FRACTIONAL_LANDUSE_IMPL_HPP

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_timing.hpp"

namespace scream {
namespace frac_landuse {

template <typename S, typename D>
std::shared_ptr<AbstractRemapper>
fracLandUseFunctions<S, D>::create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &data_file, const std::string &map_file,
    const std::string &field_name, const std::string &dim_name1,
    const std::string &dim_name2) {
  using namespace ShortFieldTagsNames;

  scorpio::register_file(data_file, scorpio::Read);
  const int ncols_data  = scorpio::get_dimlen(data_file, dim_name1);
  const int nclass_data = scorpio::get_dimlen(data_file, dim_name2);

  scorpio::release_file(data_file);

  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of
  // code)
  auto horiz_interp_tgt_grid =
      model_grid->clone("frac_land_use_horiz_interp_tgt_grid", true);

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if(ncols_data == ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(ncols_data <= ncols_model,
                     "Error! We do not allow to coarsen fractional land use "
                     "data to fit the model. We only allow\n"
                     "       fractional land use data to be at the same or "
                     "coarser resolution as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(map_file != "",
                     "ERROR: fractional land use data is on a different grid "
                     "than the model one,\n"
                     "       but spa_remap_file is missing from fractional "
                     "land use parameter list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();

  const auto tgt_grid = remapper->get_tgt_grid();

  const auto layout_2d = tgt_grid->get_2d_vector_layout(nclass_data, "class");
  const auto nondim    = ekat::units::Units::nondimensional();

  Field fractional_land_use(
      FieldIdentifier(field_name, layout_2d, nondim, tgt_grid->name()));
  fractional_land_use.allocate_view();

  remapper->register_field_from_tgt(fractional_land_use);

  remapper->registration_ends();

  return remapper;
}  // create_horiz_remapper

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
std::shared_ptr<AtmosphereInput> fracLandUseFunctions<S, D>::create_data_reader(
    const std::shared_ptr<AbstractRemapper> &horiz_remapper,
    const std::string &data_file) {
  std::vector<Field> io_fields;
  for(int i = 0; i < horiz_remapper->get_num_fields(); ++i) {
    io_fields.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(data_file, io_grid, io_fields, true);
}  // create_data_reader

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
template <typename ScalarX, typename ScalarT>
KOKKOS_INLINE_FUNCTION ScalarX fracLandUseFunctions<S, D>::linear_interp(
    const ScalarX &x0, const ScalarX &x1, const ScalarT &t) {
  return (1 - t) * x0 + t * x1;
}  // linear_interp

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
void fracLandUseFunctions<S, D>::perform_time_interpolation(
    const CommonFileRead::timeState &time_state,
    const FracLandUseInput &data_beg, const FracLandUseInput &data_end,
    const FracLandUseOutput &data_out) {
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.

  using ExeSpace = typename KT::ExeSpace;
  using ESU      = ekat::ExeSpaceUtils<ExeSpace>;

  // Gather time stamp info
  auto &t_now   = time_state.t_now;
  auto &t_beg   = time_state.t_beg_month;
  auto &delta_t = time_state.days_this_month;

  // At this stage, begin/end must have the same dimensions
  EKAT_REQUIRE(data_end.data.ncols == data_beg.data.ncols);

  const auto delta_t_fraction = (t_now - t_beg) / delta_t;

  EKAT_REQUIRE_MSG(delta_t_fraction >= 0 && delta_t_fraction <= 1,
                   "Error! Convex interpolation with coefficient out of "
                   "[0,1].\n  t_now  : " +
                       std::to_string(t_now) +
                       "\n"
                       "  t_beg  : " +
                       std::to_string(t_beg) +
                       "\n  delta_t: " + std::to_string(delta_t) + "\n");
  const int nclass  = data_beg.data.nclass;
  const int ncol    = data_beg.data.ncols;
  using ExeSpace    = typename KT::ExeSpace;
  using ESU         = ekat::ExeSpaceUtils<ExeSpace>;
  const auto policy = ESU::get_default_team_policy(ncol, nclass);

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MemberType &team) {
        const int icol = team.league_rank();  // column index
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, 0u, nclass), [&](int iclass) {
              const auto beg = data_beg.data.frac_land_use(icol, iclass);
              const auto end = data_end.data.frac_land_use(icol, iclass);
              data_out.frac_land_use(icol, iclass) =
                  linear_interp(beg, end, delta_t_fraction);
            });
      });
  Kokkos::fence();
}  // perform_time_interpolation

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
void fracLandUseFunctions<S, D>::fracLandUse_main(
    const CommonFileRead::timeState &time_state,
    const FracLandUseInput &data_beg, const FracLandUseInput &data_end,
    const FracLandUseOutput &data_out) {
  // Beg/End/Tmp month must have all sizes matching

  EKAT_REQUIRE_MSG(
      data_end.data.ncols == data_beg.data.ncols,
      "Error! FracLandUseInput data structs must have the same number of "
      "columns/levels.\n");

  // Horiz interpolation can be expensive, and does not depend on the particular
  // time of the month, so it can be done ONCE per month, *outside*
  // fracLandUse_main (when updating the beg/end states, reading them from
  // file).
  EKAT_REQUIRE_MSG(data_end.data.ncols == data_out.ncols,
                   "Error! Horizontal interpolation is performed *before* "
                   "calling fracLandUse_main,\n"
                   "       FracLandUseInput and FracLandUseOutput data structs "
                   "must have the "
                   "same number columns.\n");

  // Step 1. Perform time interpolation
  perform_time_interpolation(time_state, data_beg, data_end, data_out);
}  // fracLandUse_main

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
void fracLandUseFunctions<S, D>::update_frac_land_use_data_from_file(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    const int time_index,  // zero-based
    AbstractRemapper &horiz_interp, FracLandUseInput &input) {
  using namespace ShortFieldTagsNames;
  using ESU    = ekat::ExeSpaceUtils<typename DefaultDevice::execution_space>;
  using Member = typename KokkosTypes<DefaultDevice>::MemberType;

  start_timer("EAMxx::FracLandUse::update_frac_land_use_data_from_file");

  // 1. Read from file
  start_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::read_data");
  scorpio_reader->read_variables(time_index);
  stop_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::read_data");

  // 2. Run the horiz remapper (it is a do-nothing op if FracLandUse data is on
  // same grid as model)
  start_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::horiz_remap");
  horiz_interp.remap(/*forward = */ true);
  stop_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::horiz_remap");

  // 3. Copy from the tgt field of the remapper into the FracLandUse_data,
  // padding data if necessary
  start_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::copy_and_pad");
  // Recall, the fields are registered in the order:
  // Read the field from the file
  const auto field_from_file =
      horiz_interp.get_tgt_field(0).get_view<const Real **>();
  // copy data to the input view
  Kokkos::deep_copy(input.data.frac_land_use, field_from_file);
  Kokkos::fence();
  stop_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::copy_and_pad");

  stop_timer("EAMxx::FracLandUse::update_frac_land_use_data_from_file");

}  // END update_frac_landuse_data_from_file

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
void fracLandUseFunctions<S, D>::update_timestate(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    AbstractRemapper &horiz_interp, CommonFileRead::timeState &time_state,
    FracLandUseInput &beg, FracLandUseInput &end) {
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that FracLandUse assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1;  // Make it 0-based
  if(month != time_state.current_month) {
    // Update the FracLandUse time state information
    time_state.current_month = month;
    time_state.t_beg_month =
        util::TimeStamp({ts.get_year(), month + 1, 1}, {0, 0, 0})
            .frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(), month + 1);

    // Copy end'data into beg'data, and read in the new
    // end
    std::swap(beg, end);

    // Update the FracLandUse forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong
    // values
    //       to be assigned.  A timestep greater than a month is very unlikely
    //       so we will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    update_frac_land_use_data_from_file(scorpio_reader, ts, next_month,
                                        horiz_interp, end);
  }

}  // END updata_FracLandUse_timestate

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
void fracLandUseFunctions<S, D>::init_frac_landuse_file_read(
    const std::string field_name, const std::string dim_name1,
    const std::string dim_name2,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::string &data_file, const std::string &mapping_file,
    // output
    std::shared_ptr<AbstractRemapper> &FracLandUseHorizInterp,
    FracLandUseInput &FracLandUseData_start,
    FracLandUseInput &FracLandUseData_end,
    FracLandUseOutput &FracLandUseData_out,
    std::shared_ptr<AtmosphereInput> &FracLandUseDataReader) {
  // Init horizontal remap
  FracLandUseHorizInterp = create_horiz_remapper(
      grid, data_file, mapping_file, field_name, dim_name1, dim_name2);

  // number of columns in the file
  const int ncol = scorpio::get_dimlen(data_file, dim_name1);

  // number of fractional land use classes
  const int nclass = scorpio::get_dimlen(data_file, dim_name2);

  // Initialize the size of start/end/out data structures
  FracLandUseData_start = FracLandUseInput(ncol, nclass);
  FracLandUseData_end   = FracLandUseInput(ncol, nclass);
  FracLandUseData_out.init(ncol, nclass, true);

  // Create reader (an AtmosphereInput object)
  FracLandUseDataReader = create_data_reader(FracLandUseHorizInterp, data_file);
}  // init_frac_landuse_file_read

}  // namespace frac_landuse
}  // namespace scream

#endif  // FRACTIONAL_LANDUSE_IMPL_HPP
