#ifndef FRACTIONAL_LANDUSE_IMPL_HPP
#define FRACTIONAL_LANDUSE_IMPL_HPP

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_timing.hpp"

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
                     "       but remap file is missing from fractional "
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
void fracLandUseFunctions<S, D>::update_frac_land_use_data_from_file(
    std::shared_ptr<AtmosphereInput> &scorpio_reader,
    AbstractRemapper &horiz_interp, const_view_2d &input) {
  start_timer("EAMxx::FracLandUse::update_frac_land_use_data_from_file");

  // 1. Read from file
  start_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::read_data");
  scorpio_reader->read_variables();
  stop_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::read_data");

  // 2. Run the horiz remapper (it is a do-nothing op if FracLandUse data is on
  // same grid as model)
  start_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::horiz_remap");
  horiz_interp.remap_fwd();
  stop_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::horiz_remap");

  // 3. Get the tgt field of the remapper
  start_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::get_field");
  // Recall, the fields are registered in the order:
  // Read the field from the file
  input = horiz_interp.get_tgt_field(0).get_view<const Real **>();
  stop_timer(
      "EAMxx::FracLandUse::update_frac_land_use_data_from_file::get_field");

  stop_timer("EAMxx::FracLandUse::update_frac_land_use_data_from_file");

}  // END update_frac_landuse_data_from_file

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

template <typename S, typename D>
void fracLandUseFunctions<S, D>::init_frac_landuse_file_read(
    const int ncol, const std::string field_name, const std::string dim_name1,
    const std::string dim_name2,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::string &data_file, const std::string &mapping_file,
    // output
    std::shared_ptr<AbstractRemapper> &FracLandUseHorizInterp,
    std::shared_ptr<AtmosphereInput> &FracLandUseDataReader) {
  // Init horizontal remap
  FracLandUseHorizInterp = create_horiz_remapper(
      grid, data_file, mapping_file, field_name, dim_name1, dim_name2);

  // Create reader (an AtmosphereInput object)
  FracLandUseDataReader = create_data_reader(FracLandUseHorizInterp, data_file);
}  // init_frac_landuse_file_read

}  // namespace frac_landuse
}  // namespace scream

#endif  // FRACTIONAL_LANDUSE_IMPL_HPP
