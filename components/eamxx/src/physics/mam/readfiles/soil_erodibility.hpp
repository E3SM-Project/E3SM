#ifndef SOIL_ERODIBILITY_HPP
#define SOIL_ERODIBILITY_HPP

// For AtmosphereInput
#include "share/io/scorpio_input.hpp"

namespace scream {
namespace soil_erodibility {

template <typename ScalarType, typename DeviceType>
struct soilErodibilityFunctions {
  using Device = DeviceType;

  using KT            = KokkosTypes<Device>;
  using const_view_1d = typename KT::template view_1d<const Real>;

  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &soilErodibility_data_file, const std::string &map_file,
      const std::string &field_name, const std::string &dim_name1);

  static std::shared_ptr<AtmosphereInput> create_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &data_file);

  static void update_soil_erodibility_data_from_file(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      AbstractRemapper &horiz_interp, const_view_1d &input);

  static void init_soil_erodibility_file_read(
      const int ncol, const std::string field_name, const std::string dim_name1,
      const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file, const std::string &mapping_file,
      // output
      std::shared_ptr<AbstractRemapper> &SoilErodibilityHorizInterp,
      std::shared_ptr<AtmosphereInput> &SoilErodibilityDataReader);

};  // struct soilErodilityFunctions

}  // namespace soil_erodibility
}  // namespace scream
#endif  // SOIL_ERODIBILITY_HPP

#include "soil_erodibility_impl.hpp"