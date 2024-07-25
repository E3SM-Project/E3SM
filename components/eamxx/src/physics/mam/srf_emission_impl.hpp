#ifndef SRF_EMISSION_IMPL_HPP
#define SRF_EMISSION_IMPL_HPP

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"

namespace scream::mam_coupling {
namespace {

template <typename S, typename D>
std::shared_ptr<AbstractRemapper>
srfEmissFunctions<S, D>::create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &data_file, const std::string &map_file,
    const bool use_iop) {
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
  if(ncols_data == ncols_model or
     use_iop /*IOP class defines it's own remapper for file data*/) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(
        ncols_data <= ncols_model,
        "Error! We do not allow to coarsen srfEmiss data to fit the "
        "model. We only allow\n"
        "       srfEmiss data to be at the same or coarser resolution "
        "as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(
        map_file != "",
        "ERROR: srfEmiss data is on a different grid than the model one,\n"
        "       but srfEmiss_remap_file is missing from srfEmiss parameter "
        "list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();

  const auto tgt_grid = remapper->get_tgt_grid();

  const auto layout_2d = tgt_grid->get_2d_scalar_layout();
  const auto nondim    = ekat::units::Units::nondimensional();

  Field agr(FieldIdentifier("AGR", layout_2d, nondim, tgt_grid->name()));
  Field rco(FieldIdentifier("RCO", layout_2d, nondim, tgt_grid->name()));
  Field shp(FieldIdentifier("SHP", layout_2d, nondim, tgt_grid->name()));
  Field slv(FieldIdentifier("SLV", layout_2d, nondim, tgt_grid->name()));
  Field tra(FieldIdentifier("TRA", layout_2d, nondim, tgt_grid->name()));
  Field wst(FieldIdentifier("WST", layout_2d, nondim, tgt_grid->name()));

  agr.allocate_view();
  rco.allocate_view();
  shp.allocate_view();
  slv.allocate_view();
  tra.allocate_view();
  wst.allocate_view();

  remapper->register_field_from_tgt(agr);
  remapper->register_field_from_tgt(rco);
  remapper->register_field_from_tgt(shp);
  remapper->register_field_from_tgt(slv);
  remapper->register_field_from_tgt(tra);
  remapper->register_field_from_tgt(wst);

  remapper->registration_ends();

  return remapper;
}  // create_horiz_remapper

template <typename S, typename D>
std::shared_ptr<AtmosphereInput>
srfEmissFunctions<S, D>::create_srfEmiss_data_reader(
    const std::shared_ptr<AbstractRemapper> &horiz_remapper,
    const std::string &srfEmiss_data_file) {
  std::vector<Field> io_fields;
  for(int i = 0; i < horiz_remapper->get_num_fields(); ++i) {
    io_fields.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(srfEmiss_data_file, io_grid,
                                           io_fields, true);
}

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

  auto agr = srfEmiss_horiz_interp.get_tgt_field(0).get_view<const Real *>();
  auto rco = srfEmiss_horiz_interp.get_tgt_field(1).get_view<const Real *>();
  auto shp = srfEmiss_horiz_interp.get_tgt_field(2).get_view<const Real *>();
  auto slv = srfEmiss_horiz_interp.get_tgt_field(3).get_view<const Real *>();
  auto tra = srfEmiss_horiz_interp.get_tgt_field(4).get_view<const Real *>();
  auto wst = srfEmiss_horiz_interp.get_tgt_field(5).get_view<const Real *>();

  const auto &layout = srfEmiss_horiz_interp.get_tgt_field(0)
                           .get_header()
                           .get_identifier()
                           .get_layout();

  const int ncols = layout.dim(COL);

  auto srfEmiss_data_agr = ekat::scalarize(srfEmiss_input.data.AGR);
  auto srfEmiss_data_rco = ekat::scalarize(srfEmiss_input.data.RCO);
  auto srfEmiss_data_shp = ekat::scalarize(srfEmiss_input.data.SHP);
  auto srfEmiss_data_slv = ekat::scalarize(srfEmiss_input.data.SLV);
  auto srfEmiss_data_tra = ekat::scalarize(srfEmiss_input.data.TRA);
  auto srfEmiss_data_wst = ekat::scalarize(srfEmiss_input.data.WST);

  auto copy_and_pad = KOKKOS_LAMBDA(const Member &team) {
    int icol                = team.league_rank();
    srfEmiss_data_agr(icol) = agr(icol);
    srfEmiss_data_rco(icol) = rco(icol);
    srfEmiss_data_shp(icol) = shp(icol);
    srfEmiss_data_slv(icol) = slv(icol);
    srfEmiss_data_tra(icol) = tra(icol);
    srfEmiss_data_wst(icol) = wst(icol);
  };
  Kokkos::fence();
  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file::copy_and_pad");

  stop_timer("EAMxx::srfEmiss::update_srfEmiss_data_from_file");

}  // END update_srfEmiss_data_from_file

}  // namespace
}  // namespace scream::mam_coupling

#endif  // SRF_EMISSION_IMPL_HPP