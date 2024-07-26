#ifndef SRF_EMISSION_IMPL_HPP
#define SRF_EMISSION_IMPL_HPP

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"

namespace scream::mam_coupling {
namespace {

template <typename S, typename D>
template <std::size_t FN>
std::shared_ptr<AbstractRemapper>
srfEmissFunctions<S, D>::create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &data_file,
    const std::array<std::string, FN> &field_names,
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

  std::vector<Field> emiss_components;

  for(int icomp = 0; icomp < FN; ++icomp) {
    auto comp_name = field_names[icomp];
    // set and allocate fields
    Field f(FieldIdentifier(comp_name, layout_2d, nondim, tgt_grid->name()));
    f.allocate_view();
    emiss_components.push_back(f);
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
  std::vector<Field> emiss_components;
  for(int i = 0; i < horiz_remapper->get_num_fields(); ++i) {
    emiss_components.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(srfEmiss_data_file, io_grid,
                                           emiss_components, true);
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

  using ExeSpace = typename KT::ExeSpace;
  using ESU      = ekat::ExeSpaceUtils<ExeSpace>;

  // Gather time stamp info
  auto &t_now   = time_state.t_now;
  auto &t_beg   = time_state.t_beg_month;
  auto &delta_t = time_state.days_this_month;

  // At this stage, begin/end must have the same dimensions
  EKAT_REQUIRE(data_end.data.ncols == data_beg.data.ncols);

  auto delta_t_fraction = (t_now - t_beg) / delta_t;

  EKAT_REQUIRE_MSG(
      delta_t_fraction >= 0 && delta_t_fraction <= 1,
      "Error! Convex interpolation with coefficient out of [0,1].\n"
      "  t_now  : " +
          std::to_string(t_now) +
          "\n"
          "  t_beg  : " +
          std::to_string(t_beg) +
          "\n"
          "  delta_t: " +
          std::to_string(delta_t) + "\n");
  using KT = ekat::KokkosTypes<DefaultDevice>;
  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(
          data_beg.data.ncols, 1);

  Kokkos::parallel_for(
      "srfEmiss_time_interp_loop", policy,
      KOKKOS_LAMBDA(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) {
        const int icol = team.league_rank();

        // We have only 2d vars, so we need to make one team member handle it.
        Kokkos::single(Kokkos::PerTeam(team), [&] {
          data_out.emiss_components[1](icol) = linear_interp(
              data_beg.data.emiss_components[1](icol),
              data_end.data.emiss_components[1](icol), delta_t_fraction);

          /*data_out.AGR(icol) =
              linear_interp(data_beg.data.AGR(icol), data_end.data.AGR(icol),
                            delta_t_fraction);
          data_out.RCO(icol) =
              linear_interp(data_beg.data.RCO(icol), data_end.data.RCO(icol),
                            delta_t_fraction);
          data_out.SHP(icol) =
              linear_interp(data_beg.data.SHP(icol), data_end.data.SHP(icol),
                            delta_t_fraction);
          data_out.SLV(icol) =
              linear_interp(data_beg.data.SLV(icol), data_end.data.SLV(icol),
                            delta_t_fraction);
          data_out.TRA(icol) =
              linear_interp(data_beg.data.TRA(icol), data_end.data.TRA(icol),
                            delta_t_fraction);
          data_out.WST(icol) =
              linear_interp(data_beg.data.WST(icol), data_end.data.WST(icol),
                            delta_t_fraction);*/
        });
      });
  Kokkos::fence();
}  // perform_time_interpolation

template <typename S, typename D>
void srfEmissFunctions<S, D>::srfEmiss_main(const srfEmissTimeState &time_state,
                                            const srfEmissInput &data_beg,
                                            const srfEmissInput &data_end,
                                            const srfEmissInput &data_tmp,
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
  using ESU    = ekat::ExeSpaceUtils<typename DefaultDevice::execution_space>;
  using Member = typename KokkosTypes<DefaultDevice>::MemberType;

  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file");

  // 1. Read from file
  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::read_data");
  scorpio_reader->read_variables(time_index);
  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::read_data");

  // 2. Run the horiz remapper (it is a do-nothing op if srfEmiss data is on
  // same grid as model)
  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::horiz_remap");
  srfEmiss_horiz_interp.remap(/*forward = */ true);
  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::horiz_remap");

  // 3. Copy from the tgt field of the remapper into the srfEmiss_data, padding
  // data if necessary
  start_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::copy_and_pad");
  // Recall, the fields are registered in the order: ps, ccn3, g_sw, ssa_sw,
  // tau_sw, tau_lw

  // Get pointers for the srfEmiss_input
  /*auto srfEmiss_data_agr = ekat::scalarize(srfEmiss_input.data.AGR);
  auto srfEmiss_data_rco = ekat::scalarize(srfEmiss_input.data.RCO);
  auto srfEmiss_data_shp = ekat::scalarize(srfEmiss_input.data.SHP);
  auto srfEmiss_data_slv = ekat::scalarize(srfEmiss_input.data.SLV);
  auto srfEmiss_data_tra = ekat::scalarize(srfEmiss_input.data.TRA);
  auto srfEmiss_data_wst = ekat::scalarize(srfEmiss_input.data.WST);
*/
  const auto &layout = srfEmiss_horiz_interp.get_tgt_field(0)
                           .get_header()
                           .get_identifier()
                           .get_layout();

  const int ncols = layout.dim(COL);

  // Read fields from the file
  for(int i = 0; i < 6; ++i) {
    auto aa = srfEmiss_horiz_interp.get_tgt_field(i).get_view<const Real *>();
    Kokkos::deep_copy(srfEmiss_input.data.emiss_components[i], aa);
  }
  /*auto agr = srfEmiss_horiz_interp.get_tgt_field(0).get_view<const Real *>();
  auto rco = srfEmiss_horiz_interp.get_tgt_field(1).get_view<const Real *>();
  auto shp = srfEmiss_horiz_interp.get_tgt_field(2).get_view<const Real *>();
  auto slv = srfEmiss_horiz_interp.get_tgt_field(3).get_view<const Real *>();
  auto tra = srfEmiss_horiz_interp.get_tgt_field(4).get_view<const Real *>();
  auto wst = srfEmiss_horiz_interp.get_tgt_field(5).get_view<const Real *>();

  auto copy_and_pad = KOKKOS_LAMBDA(const Member &team) {
    int icol                = team.league_rank();
    srfEmiss_data_agr(icol) = agr(icol);
    srfEmiss_data_rco(icol) = rco(icol);
    srfEmiss_data_shp(icol) = shp(icol);
    srfEmiss_data_slv(icol) = slv(icol);
    srfEmiss_data_tra(icol) = tra(icol);
    srfEmiss_data_wst(icol) = wst(icol);
  };

  auto policy = ESU::get_default_team_policy(ncols, 1);
  Kokkos::parallel_for("", policy, copy_and_pad);*/
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
    time_state.current_month = month;
    time_state.t_beg_month =
        util::TimeStamp({ts.get_year(), month + 1, 1}, {0, 0, 0})
            .frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(), month + 1);

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

}  // namespace
}  // namespace scream::mam_coupling

#endif  // SRF_EMISSION_IMPL_HPP