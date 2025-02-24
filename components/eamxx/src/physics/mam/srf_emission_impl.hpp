#ifndef SRF_EMISSION_IMPL_HPP
#define SRF_EMISSION_IMPL_HPP

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

namespace scream::mam_coupling {
template <typename S, typename D>
std::shared_ptr<AbstractRemapper>
srfEmissFunctions<S, D>::create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &data_file, const std::vector<std::string> &sector_names,
    const std::string &map_file) {
  using namespace ShortFieldTagsNames;

  scorpio::register_file(data_file, scorpio::Read);
  const int ncols_data = scorpio::get_dimlen(data_file, "ncol");
  scorpio::release_file(data_file);

  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of
  // code)
  auto horiz_interp_tgt_grid =
      model_grid->clone("srf_emiss_horiz_interp_tgt_grid", true);

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  // if the file's grid is same as model's native grid, we identity remapper
  //  (i.e., no interpolation)
  if(ncols_data == ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(ncols_data <= ncols_model,
                     "Error! We do not allow to coarsen srfEmiss data to fit "
                     "the model. We only allow\n"
                     "srfEmiss data to be at the same or coarser resolution as "
                     "the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(
        map_file != "",
        "ERROR: srfEmiss data is on a different grid than the model one,\n"
        "but srfEmiss_remap_file is missing from srfEmiss parameter "
        "list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();

  const auto tgt_grid = remapper->get_tgt_grid();

  const auto layout_2d = tgt_grid->get_2d_scalar_layout();
  const auto nondim    = ekat::units::Units::nondimensional();

  std::vector<Field> field_emiss_sectors;

  const int sector_size = sector_names.size();
  for(int icomp = 0; icomp < sector_size; ++icomp) {
    auto comp_name = sector_names[icomp];
    // set and allocate fields
    Field f(FieldIdentifier(comp_name, layout_2d, nondim, tgt_grid->name()));
    f.allocate_view();
    field_emiss_sectors.push_back(f);
    remapper->register_field_from_tgt(f);
  }

  remapper->registration_ends();

  return remapper;
}  // create_horiz_remapper

template <typename S, typename D>
std::shared_ptr<AtmosphereInput>
srfEmissFunctions<S, D>::create_srfEmiss_data_reader(
    const std::shared_ptr<AbstractRemapper> &horiz_remapper,
    const std::string &srfEmiss_data_file) {
  std::vector<Field> field_emiss_sectors;
  for(int i = 0; i < horiz_remapper->get_num_fields(); ++i) {
    field_emiss_sectors.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(srfEmiss_data_file, io_grid,
                                           field_emiss_sectors, true);
}  // create_srfEmiss_data_reader

template <typename S, typename D>
template <typename ScalarX, typename ScalarT>
KOKKOS_INLINE_FUNCTION ScalarX srfEmissFunctions<S, D>::linear_interp(
    const ScalarX &x0, const ScalarX &x1, const ScalarT &t) {
  return (1 - t) * x0 + t * x1;
}  // linear_interp

template <typename S, typename D>
void srfEmissFunctions<S, D>::perform_time_interpolation(
    const srfEmissTimeState &time_state, const srfEmissInput &data_beg,
    const srfEmissInput &data_end, const srfEmissOutput &data_out) {
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.

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
        Real accum     = 0;
        // Parallel reduction over sectors
        // FIXME: Do we need to use Kokkos::Single for each team here???
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, nsectors),
            [&](const int i, Real &update) {
              const auto beg = data_beg.data.emiss_sectors(i, icol);
              const auto end = data_end.data.emiss_sectors(i, icol);
              update += linear_interp(beg, end, delta_t_fraction);
            },
            accum);
        // Assign the accumulated value to the output
        data_out.emiss_sectors(0, icol) = accum;
      });
  Kokkos::fence();
}  // perform_time_interpolation

template <typename S, typename D>
void srfEmissFunctions<S, D>::srfEmiss_main(const srfEmissTimeState &time_state,
                                            const srfEmissInput &data_beg,
                                            const srfEmissInput &data_end,
                                            const srfEmissOutput &data_out) {
  // Beg/End/Tmp month must have all sizes matching

  EKAT_REQUIRE_MSG(
      data_end.data.ncols == data_beg.data.ncols,
      "Error! srfEmissInput data structs must have the same number of "
      "columns/levels.\n");

  // Horiz interpolation can be expensive, and does not depend on the particular
  // time of the month, so it can be done ONCE per month, *outside*
  // srfEmiss_main (when updating the beg/end states, reading them from file).
  EKAT_REQUIRE_MSG(
      data_end.data.ncols == data_out.ncols,
      "Error! Horizontal interpolation is performed *before* "
      "calling srfEmiss_main,\n"
      "       srfEmissInput and srfEmissOutput data structs must have the "
      "same number columns.\n");

  // Step 1. Perform time interpolation
  perform_time_interpolation(time_state, data_beg, data_end, data_out);
}  // srfEmiss_main

template <typename S, typename D>
void srfEmissFunctions<S, D>::update_srfEmiss_data_from_file(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    const int time_index,  // zero-based
    AbstractRemapper &srfEmiss_horiz_interp, srfEmissInput &srfEmiss_input) {
  using namespace ShortFieldTagsNames;

  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file");

  // 1. Read from file
  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::read_data");
  scorpio_reader->read_variables(time_index);
  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::read_data");

  // 2. Run the horiz remapper (it is a do-nothing op if srfEmiss data is on
  // same grid as model)
  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::horiz_remap");
  srfEmiss_horiz_interp.remap_fwd();
  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::horiz_remap");

  // 3. Copy from the tgt field of the remapper into the srfEmiss_data, padding
  // data if necessary
  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::copy_and_pad");
  // Recall, the fields are registered in the order: ps, ccn3, g_sw, ssa_sw,
  // tau_sw, tau_lw

  // Read fields from the file
  for(int i = 0; i < srfEmiss_horiz_interp.get_num_fields(); ++i) {
    auto sector =
        srfEmiss_horiz_interp.get_tgt_field(i).get_view<const Real *>();
    const auto emiss =
        Kokkos::subview(srfEmiss_input.data.emiss_sectors, i, Kokkos::ALL());
    Kokkos::deep_copy(emiss, sector);
  }

  Kokkos::fence();
  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::copy_and_pad");

  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file");

}  // END update_srfEmiss_data_from_file

template <typename S, typename D>
void srfEmissFunctions<S, D>::update_srfEmiss_timestate(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    AbstractRemapper &srfEmiss_horiz_interp, srfEmissTimeState &time_state,
    srfEmissInput &srfEmiss_beg, srfEmissInput &srfEmiss_end) {
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that srfEmiss assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1;  // Make it 0-based
  if(month != time_state.current_month) {
    // Update the srfEmiss time state information
    time_state.current_month   = month;
    time_state.t_beg_month     = ts.curr_month_beg().frac_of_year_in_days();
    time_state.days_this_month = ts.days_in_curr_month();

    // Copy srfEmiss_end'data into srfEmiss_beg'data, and read in the new
    // srfEmiss_end
    std::swap(srfEmiss_beg, srfEmiss_end);

    // Update the srfEmiss forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong
    // values
    //       to be assigned.  A timestep greater than a month is very unlikely
    //       so we will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    update_srfEmiss_data_from_file(scorpio_reader, ts, next_month,
                                   srfEmiss_horiz_interp, srfEmiss_end);
  }

}  // END updata_srfEmiss_timestate

template <typename S, typename D>
void srfEmissFunctions<S, D>::init_srf_emiss_objects(
    const int ncol, const std::shared_ptr<const AbstractGrid> &grid,
    const std::string &data_file, const std::vector<std::string> &sectors,
    const std::string &srf_map_file,
    // output
    std::shared_ptr<AbstractRemapper> &SrfEmissHorizInterp,
    srfEmissInput &SrfEmissData_start, srfEmissInput &SrfEmissData_end,
    srfEmissOutput &SrfEmissData_out,
    std::shared_ptr<AtmosphereInput> &SrfEmissDataReader) {
  // Init horizontal remap
  SrfEmissHorizInterp =
      create_horiz_remapper(grid, data_file, sectors, srf_map_file);

  // Initialize the size of start/end/out data structures
  SrfEmissData_start = srfEmissInput(ncol, sectors.size());
  SrfEmissData_end   = srfEmissInput(ncol, sectors.size());
  SrfEmissData_out.init(ncol, 1, true);

  // Create reader (an AtmosphereInput object)
  SrfEmissDataReader =
      create_srfEmiss_data_reader(SrfEmissHorizInterp, data_file);
}  // init_srf_emiss_objects
}  // namespace scream::mam_coupling

#endif  // SRF_EMISSION_IMPL_HPP
