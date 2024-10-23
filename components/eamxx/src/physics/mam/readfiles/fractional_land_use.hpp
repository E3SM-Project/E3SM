#ifndef FRACTIONAL_LANDUSE_HPP
#define FRACTIONAL_LANDUSE_HPP

// For AtmosphereInput
#include "share/io/scorpio_input.hpp"

namespace scream {
namespace frac_landuse {

template <typename ScalarType, typename DeviceType>
struct fracLandUseFunctions {
  using Device = DeviceType;

  using KT            = KokkosTypes<Device>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  // Fractional land use routines
  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &fracLandUse_data_file, const std::string &map_file,
      const std::string &field_name, const std::string &dim_name1,
      const std::string &dim_name2);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static std::shared_ptr<AtmosphereInput> create_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &data_file);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static void update_frac_land_use_data_from_file(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      AbstractRemapper &horiz_interp, const_view_2d &input);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static void init_frac_landuse_file_read(
      const int ncol, const std::string field_name, const std::string dim_name1,
      const std::string dim_name2,
      const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file, const std::string &mapping_file,
      // output
      std::shared_ptr<AbstractRemapper> &FracLandUseHorizInterp,
      std::shared_ptr<AtmosphereInput> &FracLandUseDataReader);

};  // struct fracLandUseFunctions

}  // namespace frac_landuse
}  // namespace scream
#endif  // FRACTIONAL_LANDUSE_HPP

#include "fractional_land_use_impl.hpp"