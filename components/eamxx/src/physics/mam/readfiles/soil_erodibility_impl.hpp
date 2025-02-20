#ifndef SOIL_ERODIBILITY_IMPL_HPP
#define SOIL_ERODIBILITY_IMPL_HPP

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_timing.hpp"

namespace scream {
namespace soil_erodibility {

template <typename S, typename D>
std::shared_ptr<AbstractRemapper>
soilErodibilityFunctions<S, D>::create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &data_file, const std::string &map_file,
    const std::string &field_name, const std::string &dim_name1) {
  using namespace ShortFieldTagsNames;

  scorpio::register_file(data_file, scorpio::Read);
  const int ncols_data = scorpio::get_dimlen(data_file, dim_name1);

  scorpio::release_file(data_file);

  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of
  // code)
  auto horiz_interp_tgt_grid =
      model_grid->clone("soil_erodibility_horiz_interp_tgt_grid", true);

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if(ncols_data == ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(ncols_data <= ncols_model,
                     "Error! We do not allow to coarsen soil erodibility "
                     "data to fit the model. We only allow\n"
                     "       soil erodibility data to be at the same or "
                     "coarser resolution as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(map_file != "",
                     "ERROR: soil erodibility data is on a different grid "
                     "than the model one,\n"
                     "       but remap file is missing from soil erodibility "
                     "parameter list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();

  const auto tgt_grid = remapper->get_tgt_grid();

  const auto layout_2d = tgt_grid->get_2d_scalar_layout();
  const auto nondim    = ekat::units::Units::nondimensional();

  Field soil_erodibility(
      FieldIdentifier(field_name, layout_2d, nondim, tgt_grid->name()));
  soil_erodibility.allocate_view();

  remapper->register_field_from_tgt(soil_erodibility);

  remapper->registration_ends();

  return remapper;

}  // create_horiz_remapper

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
std::shared_ptr<AtmosphereInput>
soilErodibilityFunctions<S, D>::create_data_reader(
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
template <typename S, typename D>
void soilErodibilityFunctions<S, D>::update_soil_erodibility_data_from_file(
    std::shared_ptr<AtmosphereInput> &scorpio_reader,
    AbstractRemapper &horiz_interp, const_view_1d &input) {
  start_timer("EAMxx::soilErodibility::update_soil_erodibility_data_from_file");

  // 1. Read from file
  start_timer(
      "EAMxx::soilErodibility::update_soil_erodibility_data_from_file::read_"
      "data");
  scorpio_reader->read_variables();
  stop_timer(
      "EAMxx::soilErodibility::update_soil_erodibility_data_from_file::read_"
      "data");

  // 2. Run the horiz remapper (it is a do-nothing op if soilErodibility data is
  // on same grid as model)
  start_timer(
      "EAMxx::soilErodibility::update_soil_erodibility_data_from_file::horiz_"
      "remap");
  horiz_interp.remap_fwd();
  stop_timer(
      "EAMxx::soilErodibility::update_soil_erodibility_data_from_file::horiz_"
      "remap");

  // 3. Get the tgt field of the remapper
  start_timer(
      "EAMxx::soilErodibility::update_soil_erodibility_data_from_file::get_"
      "field");
  // Recall, the fields are registered in the order:
  // Read the field from the file
  input = horiz_interp.get_tgt_field(0).get_view<const Real *>();
  stop_timer(
      "EAMxx::soilErodibility::update_soil_erodibility_data_from_file::get_"
      "field");

  stop_timer("EAMxx::soilErodibility::update_soil_erodibility_data_from_file");

}  // END update_soil_erodibility_data_from_file

// -------------------------------------------------------------------------------------------
template <typename S, typename D>
void soilErodibilityFunctions<S, D>::init_soil_erodibility_file_read(
    const int ncol, const std::string field_name, const std::string dim_name1,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::string &data_file, const std::string &mapping_file,
    // output
    std::shared_ptr<AbstractRemapper> &soilErodibilityHorizInterp,
    std::shared_ptr<AtmosphereInput> &soilErodibilityDataReader) {
  // Init horizontal remap
  soilErodibilityHorizInterp = create_horiz_remapper(
      grid, data_file, mapping_file, field_name, dim_name1);

  // Create reader (an AtmosphereInput object)
  soilErodibilityDataReader =
      create_data_reader(soilErodibilityHorizInterp, data_file);
}  // init_soil_erodibility_file_read
}  // namespace soil_erodibility
}  // namespace scream

#endif  // SOIL_ERODIBILITY_IMPL_HPP
